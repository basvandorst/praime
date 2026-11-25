"""Blind test of the twin-class law (registration round_draft.txt).  Refuses to run if the frozen predictions or any code
file differ from the hashes below.
  usage:  test_twinclass.py            fresh decade [2e10, 2e11): reads fresh_S0.npz .. fresh_S8.npz
          test_twinclass.py --dry      seen decade [1e8, 1e9) (dry run): reads dry_S0.npz .. dry_S8.npz
Count files: tcsieve.py L0 W fresh_Ss.npz 1 with L0 = 2e10 (1 + s), W = 2e10 (dry: L0 = 1e8 (1 + s), W = 1e8).
Writes test_result.json (or test_result_dry.json)."""
import sys, os, json, math, hashlib
HERE = '/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/twinclass'
sys.path.insert(0, HERE)
DRY = '--dry' in sys.argv
HASHES = {
    'predictions_twinclass.json': '4d5f717166d256654f1d50d200be672ad580aeace93af631bdddc288c9ca1f5a',
    'predictions_dry.json': '6ce8aff99cd74d9e3ee8e7470ebd78e250e61090b7f6d63afd04719ecbfd6f55',
    'tcsieve.py': '5d00268cb849676206e00a0c8359eed9ba09762d1ba8627f9e59bacd8b54675f',
    'measure.py': '2747ca9b4044d5f467e410ae416858dac24e17a6ced0f8e2aff0d41699a59c31',
    'fourpoint.py': '7e8eb30c618cfd20e5605cb2e4fb07803eb377d7c7c08f1220ee06a4990f3e43',
}
for f, hsh in HASHES.items():
    got = hashlib.sha256(open(os.path.join(HERE, f), 'rb').read()).hexdigest()
    if got != hsh:
        sys.exit('REFUSED: %s has sha256 %s, registered %s' % (f, got, hsh))
import numpy as np
import measure as M

PRED = json.load(open(os.path.join(HERE, 'predictions_dry.json' if DRY else 'predictions_twinclass.json')))
pre = 'dry' if DRY else 'fresh'
P = PRED['predictions']; EM = PRED['error_model']; ALT = PRED['alternatives']
QS = (210, 2310, 30030)
meas = {}
tot = None
for s in range(9):
    d = np.load(os.path.join(HERE, '%s_S%d.npz' % (pre, s)))
    L0, W = PRED['windows']['S%d' % s]
    assert int(d['L0']) == L0 and int(d['W']) == W, 'count file %d does not cover the registered sub-window' % s
    arr = [d['P'].sum(0).astype(float), d['WP'].sum(0), d['X2'].sum(0).astype(float)]
    tot = arr if tot is None else [a + b for a, b in zip(tot, arr)]
    meas['S%d' % s] = {Q: M.measure(arr[0], arr[1], arr[2], Q) for Q in QS}
meas['D'] = {Q: M.measure(tot[0], tot[1], tot[2], Q) for Q in QS}

def pr(w, Q): return P['%s|%d' % (w, Q)]
def lev(w, key):
    u = {Q: meas[w][Q]['u2' + key] for Q in QS}; v = {Q: pr(w, Q)['p' + key] for Q in QS}
    return (u[210] / v[210], (134 * u[2310] - 14 * u[210]) / (134 * v[2310] - 14 * v[210]),
            (1484 * u[30030] - 134 * u[2310]) / (1484 * v[30030] - 134 * v[2310]))
res = {'mode': pre, 'windows': PRED['windows'], 'per_window': {}}
for w in ['D'] + ['S%d' % s for s in range(9)]:
    res['per_window'][w] = {'measured': {str(Q): meas[w][Q] for Q in QS},
                            'ratio_X': {str(Q): meas[w][Q]['u2X'] / pr(w, Q)['pX'] for Q in QS},
                            'ratio_R': {str(Q): meas[w][Q]['u2R'] / pr(w, Q)['pR'] for Q in QS},
                            'levels_X': lev(w, 'X'), 'levels_R': lev(w, 'R')}
subs = ['S%d' % s for s in range(9)]
# P1: class-level linear response, mod 30030
b = np.array([meas[w][30030]['beta'] for w in subs]); bh = np.array([pr(w, 30030)['beta_HL'] for w in subs])
seb = EM['beta_30030_inflation'] * math.sqrt(np.mean([pr(w, 30030)['se_beta_nominal'] ** 2 for w in subs])) / 3
P1 = {'beta_mean': float(b.mean()), 'beta_HL_mean': float(bh.mean()), 'se': seb, 'z': float((b.mean() - bh.mean()) / seb),
      'z_vs_independent': float((b.mean() - ALT['beta_independent']) / seb), 'z_vs_one_member': float((b.mean() - ALT['beta_one_member']) / seb)}
P1['pass'] = abs(P1['z']) <= 3 and abs(P1['z_vs_independent']) > 5 and abs(P1['z_vs_one_member']) > 5
# P2: fine-level remainder
l3 = np.array([res['per_window'][w]['levels_R'][2] for w in subs])
se2 = EM['R_level3_inflation'] * math.sqrt(2 / 1350) / 3
ind = []
for w in subs:
    q3, q2 = pr(w, 30030), pr(w, 2310)
    ind.append((1484 * (q3['pX'] + q3['pZ']) - 134 * (q2['pX'] + q2['pZ'])) / (1484 * q3['pR'] - 134 * q2['pR']))
P2 = {'R_level3_mean': float(l3.mean()), 'se': se2, 'z': float((l3.mean() - 1) / se2),
      'independence_prediction': float(np.mean(ind)), 'z_vs_independence': float((l3.mean() - np.mean(ind)) / se2)}
P2['pass'] = abs(P2['z']) <= 3
# P3: size of the linear part, mod 30030
rz = np.array([meas[w][30030]['u2Z'] / pr(w, 30030)['pZ'] for w in subs])
se3 = EM['u2Z_30030_inflation'] * math.sqrt(2 / 1484) / 3
pois = []
for w in subs:
    L0, W = PRED['windows'][w]; r1 = (30030 / 5760) / math.log(L0 + W / 2)
    pois.append(2 * pr(w, 30030)['kbar'] * (1 - r1) / pr(w, 30030)['pZ'])
P3 = {'u2Z_ratio_mean': float(rz.mean()), 'se': se3, 'z': float((rz.mean() - 1) / se3),
      'poisson_primes_ratio': float(np.mean(pois)), 'z_vs_poisson': float((rz.mean() - np.mean(pois)) / se3)}
P3['pass'] = abs(P3['z']) <= 3
# P4: cross part of the remainder on D (open point)
cross = {}
for Q in QS:
    ms = np.array([meas[w][Q]['mean'] for w in subs]); m = ms.sum()
    for key in ('X', 'R'):
        wm = float(np.sum(np.array([meas[w][Q]['u2' + key] for w in subs]) * ms) / m)
        wp = float(np.sum(np.array([pr(w, Q)['p' + key] for w in subs]) * ms) / m)
        cm = meas['D'][Q]['u2' + key] - wm; cp = pr('D', Q)['p' + key] - wp
        sd = 2 * wp * math.sqrt((1 - float(np.sum(ms ** 2)) / m ** 2) / (2 * (meas['D'][Q]['K'] - 1)))
        cross['%s%d' % (key, Q)] = {'cross_meas': cm, 'cross_HL': cp, 'sd': sd, 'z': (cm - cp) / sd, 'within_ratio': wm / wp}
Z = (cross['R2310']['z'] + cross['R30030']['z']) / math.sqrt(2)
P4 = {'z_2310': cross['R2310']['z'], 'z_30030': cross['R30030']['z'], 'Z': Z, 'A_twin_expected_Z': ALT['A_twin_cross_z']['joint'],
      'reading': 'CONFIRMED (excess anti-correlation)' if Z <= -3 else ('chance favoured' if Z >= -0.92 else 'OPEN')}
res.update(P1=P1, P2=P2, P3=P3, P4=P4, cross=cross,
           verdict='LAW CONFIRMED' if (P1['pass'] and P2['pass'] and P3['pass']) else 'LAW FALSIFIED')
out = 'test_result_dry.json' if DRY else 'test_result.json'
s = json.dumps(res, indent=1, default=float)
open(os.path.join(HERE, out), 'w').write(s)
print('P1 beta: %.4f vs %.4f +- %.4f (z %.2f; vs 0: %.0f se; vs 0.5: %.0f se) %s' % (P1['beta_mean'], P1['beta_HL_mean'], seb, P1['z'], P1['z_vs_independent'], P1['z_vs_one_member'], 'PASS' if P1['pass'] else 'FAIL'))
print('P2 R level 3: %.4f +- %.4f (z %.2f; independence %.3f) %s' % (P2['R_level3_mean'], se2, P2['z'], P2['independence_prediction'], 'PASS' if P2['pass'] else 'FAIL'))
print('P3 u2Z ratio 30030: %.4f +- %.4f (z %.2f; Poisson primes %.2f) %s' % (P3['u2Z_ratio_mean'], se3, P3['z'], P3['poisson_primes_ratio'], 'PASS' if P3['pass'] else 'FAIL'))
print('P4 cross part of R on D: z_2310 %.2f, z_30030 %.2f, Z %.2f -> %s' % (P4['z_2310'], P4['z_30030'], Z, P4['reading']))
for w in ['D'] + subs:
    print('  %s levels X %s  R %s' % (w, ' '.join('%.3f' % x for x in res['per_window'][w]['levels_X']), ' '.join('%.3f' % x for x in res['per_window'][w]['levels_R'])))
print('VERDICT:', res['verdict'], '  result sha256', hashlib.sha256(s.encode()).hexdigest())
