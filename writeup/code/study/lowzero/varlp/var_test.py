"""var_test.py: the registered variance test (round_draft.txt) of log|L'(rho, chi_D)| at 600 <= gamma <= 1000.

usage:  var_test.py PRED.json SET:VALDIR:ZEROFILE[,ZEROFILE] [SET:VALDIR:ZEROFILE] [--noise-only] [--boot 2000] [--band LO HI]
        SET in C, Bp; VALDIR = output of  lprime_values.py VALDIR ZEROFILE --lo 600 --hi 1000  (lfuninit to 1005);
        PRED.json = predictions_varlp.json (fresh.py).
Reuses the registered round-142 code (lprime_test.py: load, completeness rule, per-L-function statistics, design) with
its band set to [600, 1000]. Reading: OLS of the per-L-function variance of y = log|L'| - ln L on the design; the
frozen per-L-function predictions are matched through the same design (as lprime_test.py does).
"""
import sys, os, json, math, importlib.util
sys.dont_write_bytecode = True
import numpy as np
REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..'))
spec = importlib.util.spec_from_file_location('lprime_test', REPO + '/code/study/lowzero/lprime/lprime_test.py')
LT = importlib.util.module_from_spec(spec); spec.loader.exec_module(LT)
LT.LO, LT.HI = 600.0, 1000.0
if '--band' in sys.argv:                 # validation only (retrospective run on round 142's band); the registered test uses 600..1000
    i = sys.argv.index('--band'); LT.LO, LT.HI = float(sys.argv[i + 1]), float(sys.argv[i + 2])
PR = LT.PR; DEC = (2, 3, 5, 7); AMP = (2, 3, 5)
HYP = ('3pt', 'ldm', 'lead', 'none')


def readings(rows, st, idx):
    X, ns = LT.design([rows[i] for i in idx]); ip = ns + 1; H = np.linalg.pinv(X.T @ X) @ X.T
    y = np.array([st[i]['var'] for i in idx]); b = H @ y
    out = {'fit': b[ip:ip + len(PR)]}
    for h in HYP: out[h] = (H @ np.array([rows[i]['cum2_' + h] for i in idx]))[ip:ip + len(PR)]
    J = [PR.index(p) for p in AMP]
    d = out['3pt'][J] - out['lead'][J]; out['amp'] = float(((out['fit'][J] - out['lead'][J]) * d).sum() / (d * d).sum())
    d2 = out['3pt'][J] - out['ldm'][J]; out['amp_ldm'] = float(((out['fit'][J] - out['ldm'][J]) * d2).sum() / (d2 * d2).sum())
    ym = np.array([st[i]['mu'] for i in idx]); bm = H @ ym; out['mean_fit'] = bm[ip:ip + len(PR)]
    out['mean_3pt'] = (H @ np.array([rows[i]['cum1_3pt'] for i in idx]))[ip:ip + len(PR)]
    return out


def main(argv):
    noise = '--noise-only' in argv; nb = int(argv[argv.index('--boot') + 1]) if '--boot' in argv else 2000
    pred = json.load(open(argv[1])); rows, excl = [], []
    for s in [a for a in argv[2:] if ':' in a]:
        r, e = LT.load(pred, s); rows += r; excl += e
    st = [LT.stats_per_lfunction(r) for r in rows]; n = len(rows); idx = np.arange(n)
    print('L-functions kept %d; excluded %s; zeros used %d; max |L|/|L\'| %.1e' % (n, excl if excl else 'none', sum(x['n'] for x in st), max(r['chk'] for r in rows)))
    print('groups (2,3,5,7,11,13): %s' % [sum(p in r['bad'] for r in rows) for p in PR])
    R = readings(rows, st, idx); rng = np.random.default_rng(20261126)
    B = [readings(rows, st, rng.integers(0, n, n)) for _ in range(nb)]
    se = np.std([b['fit'] for b in B], axis=0); sa = float(np.std([b['amp'] for b in B])); sal = float(np.std([b['amp_ldm'] for b in B]))
    sem = np.std([b['mean_fit'] for b in B], axis=0)
    print('\nvariance reading     p: ' + ' '.join('%9d' % p for p in PR))
    print('bootstrap se            ' + ' '.join('%9.5f' % x for x in se))
    if not noise: print('fit                     ' + ' '.join('%+9.5f' % x for x in R['fit']))
    for h in HYP:
        print('pred %-5s (matched)    ' % h + ' '.join('%+9.5f' % x for x in R[h]) + ('' if noise else '   z: ' + ' '.join('%+.1f' % x for x in (R['fit'] - R[h]) / se)))
    J = [PR.index(p) for p in AMP]
    print('\npower: |pred_3pt - pred_h| / se (p = 2, 3, 5, 7):')
    for h in HYP[1:]:
        print('  vs %-5s ' % h + ' '.join('%.1f' % x for x in np.abs(R['3pt'] - R[h])[:4] / se[:4]))
    print('amplitude gamma (0 = leading order, 1 = three-point), p = 2, 3, 5: se %.4f' % sa + ('' if noise else ', reading %.4f' % R['amp']))
    print('amplitude gamma_L (0 = round-142 heuristic LDM, 1 = three-point), p = 2, 3, 5: se %.4f' % sal + ('' if noise else ', reading %.4f' % R['amp_ldm']))
    print('\nmean reading (secondary; two-point law) se ' + ' '.join('%.5f' % x for x in sem))
    print('pred 3pt (matched)     ' + ' '.join('%+9.5f' % x for x in R['mean_3pt']) + ('' if noise else '\nfit                    ' + ' '.join('%+9.5f' % x for x in R['mean_fit'])))
    if noise: return
    print('\nDECISIONS (round_draft.txt)')
    k = [PR.index(p) for p in DEC]
    ok1 = bool(np.all(np.abs(R['fit'][k] - R['3pt'][k]) <= 4 * se[k]))
    print('  V1 every p in 2,3,5,7: |fit - pred_3pt| <= 4 se: %s  (z = %s)' % (ok1, ' '.join('%+.1f' % x for x in (R['fit'][k] - R['3pt'][k]) / se[k])))
    ok2 = abs(R['amp'] - 1) <= 4 * sa
    print('  V2 amplitude gamma = %.3f +- %.3f, |gamma - 1| <= 4 se: %s' % (R['amp'], sa, ok2))
    print('  V3 leading order rejected (gamma >= 5 se and >= 0.5): %s' % (R['amp'] >= max(5 * sa, 0.5)))
    print('  V4 round-142 heuristic rejected (gamma_L >= 5 se and >= 0.5): %s  (gamma_L = %.3f +- %.3f)' % (R['amp_ldm'] >= max(5 * sal, 0.5), R['amp_ldm'], sal))
    fals = bool(np.any(np.abs(R['fit'][[0, 1]] - R['3pt'][[0, 1]]) > 6 * se[[0, 1]])) or abs(R['amp'] - 1) > 6 * sa
    print('  FALSIFIER (p = 2 or 3 off by more than 6 se, or |gamma - 1| > 6 se): %s' % fals)


if __name__ == '__main__':
    if len(sys.argv) < 3: sys.exit(__doc__)
    main(sys.argv)
