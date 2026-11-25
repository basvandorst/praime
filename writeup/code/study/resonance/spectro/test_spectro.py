"""test_spectro.py: the registered holdout test of predictions_spectro.json (refuses any other file).

usage:
  test_spectro.py sums [BUDGET_s]   compute the holdout chunk sums on (3e7, 1e9], chunks of 1e7, appended to sums_holdout.jsonl;
                                    resumable: rerun until it prints ALL DONE (each run stops after BUDGET seconds, default 540)
  test_spectro.py read              readings, decisions, alternatives -> readings_holdout.json (needs all 97 chunks)
  test_spectro.py dry [BUDGET_s]    the same code on the DISCOVERY range (0, 3e7] with the discovery windows (code check only;
                                    the frozen holdout predictions are not used), -> sums_dry.jsonl, then 'dryread'
Run every call as: (ulimit -v 1200000; timeout 600 nice .venv/bin/python -I -B test_spectro.py ...)"""
import sys, os, json, math, time, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import spectro_lib as S
import twin_lib as T

FROZEN = '7848335b687c2971a7fa3f1fa7b8b28c4d35021ce0c9e473c6412783b03adc0d'
H = S.HERE
blob = open(os.path.join(H, 'predictions_spectro.json')).read()
h = hashlib.sha256(blob.encode()).hexdigest()
if h != FROZEN:
    sys.exit('REFUSED: predictions_spectro.json has sha256 %s, registered %s' % (h, FROZEN))
PRED = json.loads(blob)
EC = ['EC:11a1', 'EC:15a1']; LAT = ['A4', 'C23']
EC_STATS = ['twinlo', 'twinhi', 'sg']; LAT_STATS = ['twinlo', 'twinhi', 'sg', 'single', 'shift2', 'shift6']
CH = 10_000_000

def config(run):
    if run == 'hold':
        V0, V1 = PRED['window_v']; L0, L1 = PRED['range']
        lines = {lab: [x['alpha'] for x in PRED['lines'][lab]['rows']] for lab in LAT}
        return dict(V=(V0, V1), VL=(V0, V1), L0=L0, L1=L1, lines=lines, fn='sums_holdout.jsonl')
    # dry: discovery range, discovery windows (twins: [sqrt 6e5, sqrt 3e7]; lines: [sqrt 1e5, sqrt 3e7])
    lines = {lab: [x['alpha'] for x in PRED['lines'][lab]['rows']] for lab in LAT}
    return dict(V=(math.sqrt(6e5), math.sqrt(3e7)), VL=(math.sqrt(1e5), math.sqrt(3e7)), L0=0, L1=30_000_000, lines=lines, fn='sums_dry.jsonl')

def run_sums(run, budget):
    cfg = config(run); fn = os.path.join(H, cfg['fn']); done = set()
    if os.path.exists(fn):
        for line in open(fn): j = json.loads(line); done.add((j['L'], j['U']))
    t0 = time.time(); left = 0
    with open(fn, 'a') as f:
        for L in range(cfg['L0'], cfg['L1'], CH):
            U = min(cfg['L1'], L + CH)
            if (L, U) in done: continue
            if time.time() - t0 > budget: left += 1; continue
            t1 = time.time()
            a = T.chunk_sums(L, U, *cfg['V'], labels=EC, stats=EC_STATS)
            if cfg['VL'] == cfg['V']:
                b = T.chunk_sums(L, U, *cfg['V'], labels=LAT, stats=LAT_STATS, lines=cfg['lines'])
            else:
                b = T.chunk_sums(L, U, *cfg['V'], labels=LAT, stats=LAT_STATS)
                c = T.chunk_sums(L, U, *cfg['VL'], labels=LAT, stats=[], lines=cfg['lines'])
                for lab in LAT: b[lab]['lines'] = c[lab]['lines']
            a.update(b)
            f.write(json.dumps({'L': L, 'U': U, 'res': a}) + '\n'); f.flush()
            print('chunk (%d, %d] %.1fs  total %.0fs' % (L, U, time.time() - t1, time.time() - t0), flush=True)
    print('ALL DONE' if left == 0 else 'budget spent: %d chunks left, run again' % left)

def load(run):
    cfg = config(run); fn = os.path.join(H, cfg['fn']); tot = {}; n = 0
    need = set((L, min(cfg['L1'], L + CH)) for L in range(cfg['L0'], cfg['L1'], CH)); have = set()
    for line in open(fn):
        j = json.loads(line)
        if (j['L'], j['U']) in have: continue
        have.add((j['L'], j['U']))
        for lab, rec in j['res'].items():
            for st, x in rec.items():
                key = (lab, st)
                if st == 'lines':
                    if key not in tot: tot[key] = dict(S=np.zeros(len(x['S']), complex), noise=0.0)
                    tot[key]['S'] += np.array([complex(*z) for z in x['S']]); tot[key]['noise'] += x['noise']
                else:
                    if key not in tot: tot[key] = dict(smooth=np.zeros(len(x['smooth']), complex), noise=0.0, count=0)
                    tot[key]['smooth'] += np.array([complex(*z) for z in x['smooth']]); tot[key]['noise'] += x['noise']; tot[key]['count'] += x['count']
    missing = need - have
    return tot, missing, cfg

def read(run):
    tot, missing, cfg = load(run)
    if missing: sys.exit('%d chunks missing; run sums first' % len(missing))
    out = {'sha256': h, 'run': run, 'stats': {}, 'lines': {}}
    V0, V1 = cfg['V']
    print('=== constellation statistics (reading = S / (C MT), expected 1 + 0i)')
    for key, P in sorted(PRED['stats'].items()):
        lab, st = P['label'], P['stat']
        MT = S.main_term_unit(lab, V0, V1) if run != 'hold' else complex(*P['MT_unit'])
        x = tot[(lab, st)]; z = x['smooth'] / (P['const'] * MT)
        sig = P['sigma_pred'] if run == 'hold' else P['discovery_sigma']
        rnd = math.sqrt(x['noise'] / 2) / abs(P['const'] * MT)
        ctrl_rms = math.sqrt(np.mean(np.abs(z[1:]) ** 2) / 2)
        out['stats'][key] = dict(reading=[z[0].real, z[0].imag], sigma_pred=sig, sigma_random=rnd, control_rms=ctrl_rms,
                                 controls=[[c.real, c.imag] for c in z[1:]], count=x['count'], alternatives=P['alternatives'])
        print('%-16s reading %+.4f %+.4fi  sigma_pred %.4f  (random %.4f, control rms %.4f)  n=%d' % (key, z[0].real, z[0].imag, sig, rnd, ctrl_rms, x['count']))
    # D1: pooled twinlo
    tw = [out['stats']['%s|twinlo' % lab] for lab in EC + LAT]
    wts = [1 / t['sigma_pred'] ** 2 for t in tw]; pooled = sum(w * t['reading'][0] for w, t in zip(wts, tw)) / sum(wts); ps = 1 / math.sqrt(sum(wts))
    D1 = abs(pooled - 1) <= 3 * ps
    # D2: models over the 12 twin readings
    keys12 = ['%s|%s' % (lab, st) for lab in EC + LAT for st in ('twinlo', 'twinhi', 'sg')]
    def chi2(model):
        s = 0.0
        for k in keys12:
            t = out['stats'][k]; e = t['alternatives'][model] if model != 'derived' else 1.0
            s += ((t['reading'][0] - e) ** 2 + t['reading'][1] ** 2) / t['sigma_pred'] ** 2
        return s
    def sep(model):
        return sum(((out['stats'][k]['alternatives'][model] - 1) / out['stats'][k]['sigma_pred']) ** 2 for k in keys12)
    c2 = {m: chi2(m) for m in ('derived', 'factorised_independence', 'no_resonance', 'integer_like_label', 'no_twin_enhancement')}
    D2 = all(c2[m] - c2['derived'] >= 25 for m in ('no_resonance', 'integer_like_label', 'no_twin_enhancement') if sep(m) >= 50)
    # D3: shifts
    sh = {k: out['stats'][k] for k in ('A4|shift2', 'A4|shift6', 'C23|shift2', 'C23|shift6')}
    c2s = {m: sum(((t['reading'][0] - (t['alternatives'][m] if m != 'derived' else 1.0)) ** 2 + t['reading'][1] ** 2) / t['sigma_pred'] ** 2 for t in sh.values())
           for m in ('derived', 'integer_like_label', 'prime_like', 'no_resonance')}
    D3 = all(abs(sh[k]['reading'][0] - 1) <= 3 * sh[k]['sigma_pred'] for k in ('A4|shift2', 'A4|shift6')) and c2s['integer_like_label'] - c2s['derived'] >= 25
    # D4: lines
    D4 = True
    for lab in LAT:
        Lp = PRED['lines'][lab]; x = tot[(lab, 'lines')]; rows = Lp['rows']
        if run == 'hold':
            Sp = np.array([complex(*r['S_pred']) for r in rows]); sig = Lp['sigma_pred']
        else:   # dry: predictions recomputed for the discovery line window
            import lines_lib as LL
            from fractions import Fraction as Fr
            u0, u1 = cfg['VL']; K1, kappa = S.K1_kappa(lab); N = S.label_info(lab)['N']
            def mts(r):
                uu = np.linspace(u0, u1, 20001); f = S.hann(uu, u0, u1) * uu ** 0.5 * (1 - 1j * kappa * math.sqrt(N) / (S.TWO_PI * uu * math.sqrt(r)))
                hh = uu[1] - uu[0]; s = f[0] + f[-1] + 4 * f[1:-1:2].sum() + 2 * f[2:-1:2].sum(); return 1.5 * K1 * r ** -0.25 * s * hh / 3
            Sp = np.array([r['c'] * mts(float(Fr(r['r']))) for r in rows]); sig = math.sqrt(x['noise'] / 2)
        Sm = x['S']; top = np.array([r['role'] == 'top60' for r in rows])
        slope = float(np.real(np.sum(Sm[top] * np.conj(Sp[top]))) / np.sum(np.abs(Sp[top]) ** 2)); se = sig / math.sqrt(float(np.sum(np.abs(Sp[top]) ** 2)))
        dev = np.abs(Sm - Sp) / sig
        ok_slope = abs(slope - 1) <= max(3 * se, 0.02)
        ok_extra = all(dev[i] <= 4 for i, r in enumerate(rows) if r['role'] == 'extra')
        ok_zero = all(abs(Sm[i]) / sig <= 4 for i, r in enumerate(rows) if r['role'] == 'zero')
        D4 = D4 and ok_slope and ok_extra and ok_zero
        chi2_lines = float(np.sum(dev[top] ** 2)); chi2_null = float(np.sum((np.abs(Sm[top]) / sig) ** 2))
        out['lines'][lab] = dict(slope=slope, slope_sigma=se, sigma=sig, sigma_random=math.sqrt(x['noise'] / 2), chi2_top60=chi2_lines, chi2_no_lines=chi2_null,
                                 rows=[dict(r=r['r'], role=r['role'], S_meas=[Sm[i].real, Sm[i].imag], S_pred=[Sp[i].real, Sp[i].imag], dev=float(dev[i])) for i, r in enumerate(rows)])
        print('lines %s: slope %.4f +- %.4f  chi2(top 60) %.1f / 120 (no lines %.1f); extras %s; zeros |S|/sig %s' % (lab, slope, se, chi2_lines, chi2_null,
              ['%s dev %.2f' % (r['r'], dev[i]) for i, r in enumerate(rows) if r['role'] == 'extra'],
              ['%s %.2f' % (r['r'], abs(Sm[i]) / sig) for i, r in enumerate(rows) if r['role'] == 'zero']))
    fals = not (0.8 <= pooled <= 1.2) or not (0.9 <= out['stats']['A4|shift2']['reading'][0] <= 1.1) or \
        any(not (0.95 <= out['lines'][lab]['slope'] <= 1.05) for lab in LAT)
    out['decisions'] = dict(D1=dict(pooled_twinlo=float(pooled), sigma=float(ps), passed=bool(D1)),
                            D2=dict(chi2={m: float(v) for m, v in c2.items()}, separation={m: float(sep(m)) for m in c2 if m != 'derived'}, passed=bool(D2)),
                            D3=dict(chi2={m: float(v) for m, v in c2s.items()}, passed=bool(D3)), D4=dict(passed=bool(D4)), falsified=bool(fals))
    print('D1 pooled twinlo %.4f +- %.4f  -> %s' % (pooled, ps, 'PASS' if D1 else 'FAIL'))
    print('D2 chi2 over 12 twin readings', {m: round(v, 1) for m, v in c2.items()}, 'separations', {m: round(sep(m), 1) for m in c2 if m != 'derived'}, '->', 'PASS' if D2 else 'FAIL')
    print('D3 shift chi2', {m: round(v, 1) for m, v in c2s.items()}, '->', 'PASS' if D3 else 'FAIL')
    print('D4 lines ->', 'PASS' if D4 else 'FAIL')
    print('FALSIFIED' if fals else 'not falsified')
    if run != 'hold': print('(dry run: sigma_pred replaced by the discovery sigma; decisions printed for the code check only)')
    json.dump(out, open(os.path.join(H, 'readings_%s.json' % ('holdout' if run == 'hold' else 'dry')), 'w'), indent=1)

if __name__ == '__main__':
    mode = sys.argv[1]; budget = float(sys.argv[2]) if len(sys.argv) > 2 else 540.0
    if mode == 'sums': run_sums('hold', budget)
    elif mode == 'read': read('hold')
    elif mode == 'dry': run_sums('dry', budget)
    elif mode == 'dryread': read('dry')
    else: sys.exit(__doc__)
