"""CM pilot (retrospective, 6 curves, not registered): measured Sigma2 against the law (deg2cm.CurveAny: additive primes
are missing orbits, Sym^2 bad factors from PARI) and against the alternatives.
usage: cm_pilot.py SET.json ZERODIR"""
import sys, os, json, math, pickle
import numpy as np, cypari2
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import deg2lib as D, deg2cm as CM
pari = cypari2.Pari(); pari.allocatemem(400 * 10 ** 6)
H = os.path.dirname(os.path.abspath(__file__)); Z = pickle.load(open(os.path.join(H, 'zgrid.pkl'), 'rb'))
S = json.load(open(sys.argv[1])); lo, hi = S['band']; LS = [0.5 * i for i in range(1, 21)]; IX = [int(round(l / D.W)) - 1 for l in LS]
out = []
def rep(s): print(s, flush=True); out.append(s)
show = [LS.index(l) for l in (1.0, 2.0, 4.0, 6.0, 8.0, 10.0)]
rep('l shown: ' + ' '.join('%7.1f' % LS[i] for i in show))
R = {k: [] for k in ('meas', 'full', 'k1', 'deg1', 'gue')}
for c in S['curves']:
    cv = CM.CurveAny(pari, c['ainv'])
    g = np.sort(np.load(os.path.join(sys.argv[2], c['id'] + '.npy'))); b = g[(g >= lo) & (g <= hi)]
    u = D.theta_E(c['N'], b) / math.pi; meas = []
    for l in LS:
        x = np.arange(u[0], u[-1] - l, 0.01); cnt = np.searchsorted(u, x + l) - np.searchsorted(u, x); meas.append(np.var(cnt))
    meas = np.array(meas)
    full = D.sigma2_curve(Z, cv, lo, hi)[IX]; k1 = D.sigma2_curve(Z, cv, lo, hi, kmode='k1')[IX]
    deg1 = D.sigma2_deg1(Z, cv.N, cv.bad, lo, hi)[IX]; gue = D.sigma2_gue()[IX]
    for k, v in (('meas', meas), ('full', full), ('k1', k1), ('deg1', deg1), ('gue', gue)): R[k].append(v)
    inert = sum(1 for p in D.primes_upto(50) if int(p) not in cv.bad and cv.ap[int(p)] == 0)
    rep('%s N=%d j=%d additive %s, Sym2 bad factors %s; a_p = 0 for %d of the good p < 50' % (c['id'], cv.N, c['j'], cv.add,
        {p: [round(x, 4) for x in v] for p, v in cv.sympoly.items()}, inert))
    for k, v in (('measured', meas), ('law', full), ('k>=2 flat', k1), ('deg1', deg1), ('GUE', gue)):
        rep('   %-10s ' % k + ' '.join('%7.4f' % v[i] for i in show))
M = {k: np.array(v) for k, v in R.items()}
rep('mean |measured - law| %.4f, mean |measured - k1| %.4f, mean |measured - deg1| %.4f (all 20 lengths, 6 curves)' % (
    np.abs(M['meas'] - M['full']).mean(), np.abs(M['meas'] - M['k1']).mean(), np.abs(M['meas'] - M['deg1']).mean()))
for alt in ('k1', 'deg1'):
    x = (M['full'] - M[alt]).ravel(); y = (M['meas'] - M[alt]).ravel()
    rep('slope of (measured - %s) on (law - %s), uncentred, all cells: %.3f' % (alt, alt, float((x * y).sum() / (x * x).sum())))
open(os.path.join(H, 'cm_pilot_out.txt'), 'w').write('\n'.join(out) + '\n')
