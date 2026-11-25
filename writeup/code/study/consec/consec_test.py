"""consec_test.py: registered decisions of the consecutive-primes test.  Refuses any predictions file but the frozen one.

usage: consec_test.py hold READINGS_MAIN.json READINGS_EC.json [OUT.json]     (registered run on the holdout)
       consec_test.py dry  READINGS_DISC.json                                (dry run on the discovery range; not blind)
predictions_consec.json must sit next to this script and have sha256 = FROZEN."""
import sys, os, json, math, hashlib
from scipy.stats import chi2 as CHI2

HERE = os.path.dirname(os.path.abspath(__file__))
FROZEN = 'f3e3599d43a6d5d9e73eefa1a46a1a96686ca88c58b3ce07623188096b4f3cb1'
PRED = os.path.join(HERE, 'predictions_consec.json')
h = hashlib.sha256(open(PRED, 'rb').read()).hexdigest()
if h != FROZEN: sys.exit('REFUSED: predictions_consec.json has sha256 %s, not the frozen %s' % (h, FROZEN))
P = json.load(open(PRED))
mode = sys.argv[1]
if mode == 'hold':
    RM, RE = json.load(open(sys.argv[2])), json.load(open(sys.argv[3])); outf = sys.argv[4] if len(sys.argv) > 4 else 'test_hold_out.json'
    for run, Rd in (('main', RM), ('ec', RE)):
        want = P['runs'][run]['range']
        if Rd['range'] != want: sys.exit('REFUSED: readings for run %s cover %s, registered %s' % (run, Rd['range'], want))
        if Rd['nprimes'] != P['runs'][run]['nprimes_expected']: print('WARNING: %s primes %d, expected %d' % (run, Rd['nprimes'], P['runs'][run]['nprimes_expected']))
    SRC = {nm: RM for nm in P['runs']['main']['schemes']}; SRC.update({nm: RE for nm in P['runs']['ec']['schemes']})
elif mode == 'dry':
    RD = json.load(open(sys.argv[2])); outf = 'test_dry_out.json'
    SRC = {nm: RD for run in P['runs'].values() for nm in run['schemes']}
else: sys.exit(__doc__)
res = {'mode': mode, 'sha256': h, 'decisions': {}, 'secondary': {}, 'watch': {}}
def dec(name, ok, detail, falsified=False):
    res['decisions'][name] = dict(passed=bool(ok), falsifies_H0=bool(falsified), detail=detail)
    print('%-4s %s%s  %s' % (name, 'PASS' if ok else 'FAIL', '  ** H0 FALSIFIED **' if falsified else '', detail))
names = list(SRC)
# D1a
tot = sum(SRC[nm][nm]['indep'][k]['chi2'] for nm in names for k in ('1', '2', '3'))
dof = sum(SRC[nm][nm]['indep'][k]['dof'] for nm in names for k in ('1', '2', '3')); p = float(CHI2.sf(tot, dof))
dec('D1a', p >= 1e-3, 'pooled contingency chi2 %.1f on %d dof, p = %.4g' % (tot, dof, p), falsified=p < 1e-5)
# D1b
cells = [(nm, k, v['z']) for nm in names for k, v in SRC[nm][nm]['corr'].items()]
mx = max(cells, key=lambda c: abs(c[2]))
dec('D1b', abs(mx[2]) < 4.5, '%d cells, sum z^2 = %.1f, max |z| = %.2f (%s %s)' % (len(cells), sum(c[2] ** 2 for c in cells), abs(mx[2]), mx[0], mx[1]), falsified=abs(mx[2]) >= 5)
# D1c
ws = [(nm, SRC[nm][nm]['w4']['within'], SRC[nm][nm]['w4']['within_err']) for nm in names]
wsum = sum(1 / e ** 2 for _, _, e in ws); wm = sum(v / e ** 2 for _, v, e in ws) / wsum; we = 1 / math.sqrt(wsum)
each = max(abs(v - 1) / e for _, v, e in ws)
dec('D1c', abs(wm - 1) / we < 4 and each < 4.5, 'pooled within %.4f +- %.4f (z %+.2f); max per-scheme |z| %.2f; ' % (wm, we, (wm - 1) / we, each)
    + ' '.join('%s %.4f+-%.4f' % w for w in ws), falsified=abs(wm - 1) / we > 5)
# D2
bad = []
for nm in names:
    for f, g in SRC[nm][nm]['gap'].items():
        if g['p'] < 1e-3: bad.append('%s %s p=%.2g' % (nm, f, g['p']))
        if 'gap_mult' in g and g['gap_mult']['z'] == g['gap_mult']['z'] and abs(g['gap_mult']['z']) >= 4: bad.append('%s %s gap=0 mod q z=%.2f' % (nm, f, g['gap_mult']['z']))
dec('D2', not bad, '; '.join(bad) if bad else 'all gap chi2 p >= 1e-3, all gap = 0 mod q contrasts |z| < 4')
# D3
nb = [(nm, k, v['z']) for nm in names for k, v in SRC[nm][nm]['nbr'].items()]
mx = max(nb, key=lambda c: abs(c[2])) if nb else ('-', '-', 0.0)
dec('D3', abs(mx[2]) < 4.5, '%d cells, max |z| = %.2f (%s %s)' % (len(nb), abs(mx[2]), mx[0], mx[1]))
# D4, D5
r4 = [(nm, SRC[nm][nm]['runs']['total']['p']) for nm in names if 'runs' in SRC[nm][nm]]
dec('D4', min(p for _, p in r4) >= 1e-3, ' '.join('%s p=%.3g' % x for x in r4))
r5 = [(nm, SRC[nm][nm]['triple']['p']) for nm in names if 'triple' in SRC[nm][nm]]
dec('D5', min(p for _, p in r5) >= 1e-3, ' '.join('%s p=%.3g' % x for x in r5))
# D6: short-shift dispersion
bad, pooled = [], {}
for nm in names:
    R0, R1 = SRC[nm]['range']
    for H in (100, 1000):
        ent = SRC[nm][nm]['disp']['H%d_%d_%d' % (H, R0, R1)]
        for f, v in ent.items():
            if f == 'chi' or v['R'] != v['R']: continue
            z = (v['R'] - 1) / v['err']; pooled.setdefault(H, []).append((v['R'], v['err']))
            if abs(z) >= 4: bad.append('%s H%d %s R=%.4f z=%.1f' % (nm, H, f, v['R'], z))
txt = []
for H, l in pooled.items():
    ws_ = sum(1 / e ** 2 for _, e in l); m = sum(r / e ** 2 for r, e in l) / ws_; e = 1 / math.sqrt(ws_); txt.append('H=%d pooled R %.5f +- %.5f (z %+.2f)' % (H, m, e, (m - 1) / e))
    if abs(m - 1) / e >= 4: bad.append('pooled H%d z=%.1f' % (H, (m - 1) / e))
dec('D6', not bad, '; '.join(txt + bad))
falsified = any(d['falsifies_H0'] for d in res['decisions'].values())
allpass = all(d['passed'] for d in res['decisions'].values())
res['H0'] = 'FALSIFIED' if falsified else ('STANDS' if allpass else 'STANDS (open anomaly in a non-falsifying decision)')
print('H0:', res['H0'])
# ------------------------------------------------------------------ secondary (known laws; not part of the verdict)
for key, s in P['secondary'].items():
    if '|abelian' not in key: continue
    nm = key.split('|')[0]; E = SRC[nm][nm]
    if mode == 'dry': continue                                   # LOS numbers are frozen for the holdout ranges only
    rat = {k: E['ab_b'][k]['b'] / s['b_LOS'][k] for k in ('1', '2', '3')}
    w = E.get('w4a', {'within': float('nan'), 'within_err': float('nan')})
    okw = (abs(w['within'] - s['w4a_within_LOS']) < 4 * w['within_err'] + 0.015) if 'w4a' in E else None
    res['secondary'][key] = dict(b_ratio=rat, b_ok=all(0.6 <= r <= 1.5 for r in rat.values()), w4a=[w['within'], w['within_err'], s['w4a_within_LOS']], w4a_ok=okw)
    print('  %s: b measured/LOS %s  abelian within %.4f+-%.4f (LOS %.4f) %s' % (key, ' '.join('%.2f' % r for r in rat.values()), w['within'], w['within_err'], s['w4a_within_LOS'], {True: 'ok', False: 'MISS', None: 'n/a'}[okw]))
if mode == 'hold':
    for run in ('main', 'ec'):
        bks = P['secondary']['bks_' + run]; nmiss, ntot = 0, 0
        for key, s in bks.items():
            nm, hk, f = key.split('|'); H = int(hk.split('_')[0][1:])
            if H not in (100000, 1000000): continue
            v = SRC[nm][nm]['disp'].get(hk, {}).get(f)
            if v is None or v['blocks'] < 50 or v['R'] != v['R']: continue
            ok = abs(v['R'] - s['law']) < 4 * math.sqrt(v['err'] ** 2 + 0.03 ** 2); ntot += 1; nmiss += not ok
            res['secondary'][key] = dict(R=v['R'], err=v['err'], law=s['law'], ok=ok)
        print('  BKS column law (%s, H = 1e5, 1e6, bins with >= 50 blocks): %d of %d cells within the band' % (run, ntot - nmiss, ntot))
    for w, s in P['watch'].items():
        if w == 'W1': v = SRC['C31']['C31']['runs']['1']; val = (v['chi2'], v['dof'], float(CHI2.sf(v['chi2'], v['dof'])))
        elif w == 'W2': v = SRC['C23']['C23']['triple']; val = (v['chi2'], v['dof'], v['p'])
        elif w == 'W3': val = (SRC['E37m3']['E37m3']['corr']['k2_f0_f0']['z'],)
        elif w == 'W4': v = SRC['G']['G']['runs']['1']; val = (v['chi2'], v['dof'], float(CHI2.sf(v['chi2'], v['dof'])))
        res['watch'][w] = val; print('  watch %s: %s  (%s)' % (w, val, s['rule']))
json.dump(res, open(os.path.join(HERE, outf), 'w'), indent=1)
