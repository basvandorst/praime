"""Exploratory (post hoc, not registered): what form does the genus-2 excess take? Least squares on the real residual
Re(G - pred) per curve on (10, 150] (the lowest zeros carry a separate offset), with jackknife errors over curves:
  M1 a (additive), M2 delta * (-pred) (proportional), M3 a + delta * (-pred), M4 a + a1 eps + delta * (-pred).
Also: residual variance against c, and the alternative constants c * L5^k (L5 = partial Euler product of L(1, Lambda^2_0)
over good p <= 1000, from exact local factors), for the 51 curves with exact local data."""
import sys, os, math, glob; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from g2data import load, slope
rows = load(); a0, b0 = 10, 150
for r in rows:
    m = (r['gg'] > a0) & (r['gg'] <= b0); r['Gw'] = complex(r['ph'][m].sum()); r['pw'] = -r['eps'] * r['c'] * (b0 - a0) / (2 * math.pi)
    r['res'] = (r['Gw'] - r['pw']).real
y = np.array([r['res'] for r in rows]); P = np.array([r['pw'] for r in rows]); E = np.array([r['eps'] for r in rows])
designs = {'M1 a': np.vstack([np.ones_like(y)]).T, 'M2 delta': np.vstack([-P]).T, 'M3 a + delta': np.vstack([np.ones_like(y), -P]).T,
           'M4 a + a1 eps + delta': np.vstack([np.ones_like(y), E, -P]).T}
n = len(y)
for nm, X in designs.items():
    beta = np.linalg.lstsq(X, y, rcond=None)[0]; rss = float(((y - X @ beta) ** 2).sum())
    loo = np.array([np.linalg.lstsq(X[np.arange(n) != i], y[np.arange(n) != i], rcond=None)[0] for i in range(n)])
    se = np.sqrt((n - 1) / n * ((loo - loo.mean(0)) ** 2).sum(0)); k = X.shape[1]
    print('%-22s coef %s  se %s  rss %.0f  AIC %.1f' % (nm, np.array2string(beta, precision=3), np.array2string(se, precision=3), rss, n * math.log(rss / n) + 2 * k))
cc = np.array([r['c'] for r in rows])
for lo, hi in ((-1, 1), (1, 2.5), (2.5, 6)):
    m = (cc > lo) & (cc <= hi); print('residual sd on (10,150] for c in (%g, %g]: %.1f (n = %d)' % (lo, hi, y[m].std(), m.sum()))
# alternative constants with L5 from exact factors p <= 1000
L5 = {}
for fn in glob.glob('cown/*.npz'):
    z = np.load(fn); key = os.path.basename(fn)[:-4]; p = z['p']; N = int(z['N']); e1, e2 = z['e1'], z['e2']
    lg = 0.0
    for q, a1, a2 in zip(p[p <= 1000], e1[p <= 1000], e2[p <= 1000]):
        if N % int(q) == 0: continue
        s = math.sqrt(max(a1 * a1 - 4 * (a2 - 2), 0.0)); t1, t2 = (a1 + s) / 2, (a1 - s) / 2      # 2 cos th_i
        c1, c2 = t1 / 2, t2 / 2; c1 = max(-1, min(1, c1)); c2 = max(-1, min(1, c2))
        th1, th2 = math.acos(c1), math.acos(c2); X = 1 / q
        det = (1 - X) * (1 - 2 * math.cos(th1 + th2) * X + X * X) * (1 - 2 * math.cos(th1 - th2) * X + X * X)
        lg -= math.log(det)
    L5[key] = math.exp(lg)
S = [r for r in rows if r['key'] in L5]
v = np.array([L5[r['key']] for r in S]); print('\nL5 partial (p <= 1000) over %d curves: range %.3f..%.3f, mean %.3f' % (len(S), v.min(), v.max(), v.mean()))
for k in (-1, -0.5, 0, 0.5, 1):
    Gs = [r['G'] for r in S]; Pk = [-r['eps'] * r['c'] * L5[r['key']] ** k * 150 / (2 * math.pi) for r in S]
    b, se = slope(Gs, Pk); bb = b
    rss = sum((g.real - bb * p) ** 2 for g, p in zip(Gs, Pk))
    print('constant c * L5^%+.1f : slope %.3f +- %.3f, rss with fitted slope %.0f' % (k, b, se, rss))
