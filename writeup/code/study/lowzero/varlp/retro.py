"""retro.py: RETROSPECTIVE comparison of the three-point prediction with round 142 (sets R120, Dn, Fn; heights 100-600;
D = 1153 excluded as there). For every L-function: the zero-weighted average over the band of dVar(S, L(t)) and of
dmu(S, L(t)) (assemble.Theory), then matched to the registered design (lprime_test.design: intercept, set dummies,
mean ln L, indicators of 2, 3, 5, 7, 11, 13) exactly as lprime_test.py matches its frozen predictions.
Also, retrospectively, the per-set fits from round 142's per-zero values (writeup/data/study/r142/vals_*.npz).
usage: retro.py OUT.json
"""
import sys, os, json, math, importlib.util
sys.dont_write_bytecode = True
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from assemble import Theory
REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..'))


def load_mod(name, path):
    spec = importlib.util.spec_from_file_location(name, path); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


LT = load_mod('lprime_test', REPO + '/code/study/lowzero/lprime/lprime_test.py')
LP = load_mod('lprime_predict', REPO + '/code/study/lowzero/lprime/lprime_predict.py')
PR = LT.PR
MEAS = {'var': ([-0.19422, -0.10706, -0.05408, -0.03256, -0.01385, -0.00886], [0.00160, 0.00172, 0.00197, 0.00237, 0.00344, 0.00340]),
        'mean': ([0.04327, 0.03885, 0.03007, 0.02356, 0.01636, 0.01346], [0.00073, 0.00073, 0.00091, 0.00092, 0.00158, 0.00147])}
REVIEW_MEAN = [0.0428, 0.0374, 0.0296, 0.0236]
X_, W_ = np.polynomial.legendre.leggauss(24)


def band_avg(th, D, bad, lo, hi):
    q = abs(D); a = 0 if D > 0 else 1
    t = lo + (hi - lo) * (X_ + 1) / 2; w = (hi - lo) / 2 * W_ * LP.density(t, q, a); w = w / w.sum()
    Lt = np.log(q * t / (2 * math.pi)); Lc = np.clip(Lt, 7.0, 12.5)
    dv = np.array([th.dvar(bad, L) if bad else 0.0 for L in Lc]); dm = np.array([th.dvar(bad, L, parts=True)['dmu'] if bad else 0.0 for L in Lc])
    return float((w * dv).sum()), float((w * dm).sum()), float(Lt.min()), float(Lt.max())


def main(out):
    th = Theory(); P = json.load(open(REPO + '/code/study/lowzero/lprime/predictions_lprime.json'))
    rows = []
    for name in ('R120', 'Dn', 'Fn'):
        for e in P[name]:
            if e['D'] == 1153: continue
            dv, dm, l0, l1 = band_avg(th, e['D'], e['bad'], 100.0, 600.0)
            rows.append(dict(e, set=name, cum2_3pt=dv, cum1_3pt=dm, Lmin=l0, Lmax=l1))
    print('L-functions %d; L range %.2f .. %.2f' % (len(rows), min(r['Lmin'] for r in rows), max(r['Lmax'] for r in rows)))
    X, ns = LT.design(rows); ip = ns + 1; H = np.linalg.pinv(X.T @ X) @ X.T
    match = lambda key: (H @ np.array([r[key] for r in rows]))[ip:ip + len(PR)]
    res = {'n': len(rows), 'primes': list(PR)}
    for stat, k3, kl in (('var', 'cum2_3pt', 'cum2'), ('mean', 'cum1_3pt', 'cum1')):
        m3 = match(k3); fit, se = map(np.array, MEAS[stat]); res[stat] = {'three_point': list(m3), 'measured': list(fit), 'se': list(se)}
        print('\n%s (RETROSPECTIVE; round 142 fit and bootstrap se from test_F3_out.txt)' % stat)
        print('  p:                ' + ' '.join('%9d' % p for p in PR))
        print('  measured          ' + ' '.join('%+9.5f' % x for x in fit))
        print('  three-point       ' + ' '.join('%+9.5f' % x for x in m3) + '   z: ' + ' '.join('%+.1f' % x for x in (fit - m3) / se))
        for h in ('ldm', 'lead'):
            mh = match(kl + '_' + h); res[stat][h] = list(mh)
            print('  %-16s  ' % h + ' '.join('%+9.5f' % x for x in mh) + '   z: ' + ' '.join('%+.1f' % x for x in (fit - mh) / se))
        if stat == 'var':
            J = [0, 1, 2]; lead = np.array(res['var']['lead'])[J]; w = 1 / se[J] ** 2
            g = float((w * (fit[J] - lead) * (m3[J] - lead)).sum() / (w * (m3[J] - lead) ** 2).sum())
            res['var']['amplitude_vs_lead'] = g
            print('  amplitude of the three-point 1/L-and-beyond term (fit - lead)/(3pt - lead), p = 2, 3, 5: %.3f' % g)
        else:
            print('  review48 (swapped term incl.) ' + ' '.join('%+.4f' % x for x in REVIEW_MEAN))
    # per-set (retrospective, from the stored per-zero values)
    print('\nper set (RETROSPECTIVE: per-L-function variance of log|L\'| - ln L from r142/vals_*.npz, OLS on [1, mean ln L, indicators])')
    res['per_set'] = {}
    for name in ('R120', 'Dn', 'Fn'):
        z = np.load(REPO + '/data/study/r142/vals_%s.npz' % name); R = [r for r in rows if r['set'] == name]; Y = []; keep = []
        for r in R:
            k = 'D%d' % r['D']
            if k + '__gamma' not in z.files: continue
            g = z[k + '__gamma']; y = z[k + '__logabs'] - np.log(np.log(abs(r['D']) * g / (2 * math.pi))); sel = (g >= 100) & (g <= 600)
            Y.append(y[sel].var(ddof=1)); keep.append(r)
        Xs = np.array([[1.0, r['mean_lnL']] + [1.0 * (p in r['bad']) for p in PR] for r in keep]); Xs[:, 1:] -= Xs[:, 1:].mean(0)
        Hs = np.linalg.pinv(Xs.T @ Xs) @ Xs.T; b = Hs @ np.array(Y); m3 = Hs @ np.array([r['cum2_3pt'] for r in keep])
        ml = Hs @ np.array([r['cum2_ldm'] for r in keep]); res['per_set'][name] = dict(n=len(keep), fit=list(b[2:]), three_point=list(m3[2:]), ldm=list(ml[2:]),
                                                                                        meanL=float(np.mean([math.exp(r['mean_lnL']) for r in keep])))
        print('  %-4s n=%3d  <L>=%.2f  fit  %s' % (name, len(keep), res['per_set'][name]['meanL'], ' '.join('%+.4f' % x for x in b[2:6])))
        print('                       3pt  %s   ldm %s' % (' '.join('%+.4f' % x for x in m3[2:6]), ' '.join('%+.4f' % x for x in ml[2:6])))
    json.dump(res, open(out, 'w'), indent=1)
    keep = ('D', 'bad', 'mean_lnL', 'n_zeros_expected', 'cum2_3pt', 'cum1_3pt', 'cum2_ldm', 'cum2_lead', 'cum2_none')
    json.dump({s: [{k: r[k] for k in keep} for r in rows if r['set'] == s] for s in ('R120', 'Dn', 'Fn')}, open('retro_pred_r142.json', 'w'))


if __name__ == '__main__':
    main(sys.argv[1])
