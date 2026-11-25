import json, math, sys, glob
import numpy as np
sys.path.insert(0, '.')
from r1_offsets import theta, jk_mean, load
def res(g, Q, mu, eps, c):
    gg = g[(g > 0) & (g <= 150)]; return np.exp(2j * theta(gg, Q, mu)).sum().real + eps * c * 150 / (2 * math.pi)
R = []
for fn in glob.glob('rep/*.npz'):
    z = np.load(fn); i = json.loads(str(z['info'])); R.append((int(i['key'].split('_')[-1]), res(z['g'], i['Q'], i['mu'], i['eps'][0], i['c']), i['eps'][0], i['c']))
R.sort(); print('idx 0-5 residuals', [round(r[1], 2) for r in R[:6]])
old = [r for r in load() if r['cls'] == 'g2']
O = [(res(r['g'], r['Q'], r['mu'], r['eps'], r['c']), r['eps'], r['c']) for r in old]
U = 150 / (2 * math.pi)
def fits(rows, lab):
    x = np.array([r[0] for r in rows]); e = np.array([r[1] for r in rows]); c = np.array([r[2] for r in rows]); n = len(x)
    for nm, X in (('a', np.c_[np.ones(n)]), ('b*eps c U', np.c_[e * c * U]), ('a + b*eps c U', np.c_[np.ones(n), e * c * U])):
        b, *_ = np.linalg.lstsq(X, x, rcond=None); r = x - X @ b
        # jackknife se of coefficients
        loo = np.array([np.linalg.lstsq(np.delete(X, k, 0), np.delete(x, k), rcond=None)[0] for k in range(n)])
        se = np.sqrt((n - 1) / n * ((loo - loo.mean(0)) ** 2).sum(0))
        print('  %-15s %-14s coef %s se %s  rss %.0f' % (lab, nm, np.round(b, 3), np.round(se, 3), r @ r))
fits(O, 'rounds 151+152'); fits([(r[1], r[2], r[3]) for r in R], 'replication'); fits(O + [(r[1], r[2], r[3]) for r in R], 'pooled 100')
# eps+1 slope with offset 1.41 removed
for lab, rows in (('replication', [(r[1], r[2], r[3]) for r in R]), ('pooled 100', O + [(r[1], r[2], r[3]) for r in R])):
    x = np.array([r[0] for r in rows]); e = np.array([r[1] for r in rows]); c = np.array([r[2] for r in rows]); P = -e * c * U
    for off in (0.0, 1.41):
        G = x - off + P
        for nm, sel in (('all', e != 0), ('eps+1', e > 0), ('eps-1', e < 0)):
            a = G[sel] * P[sel]; b = P[sel] ** 2; sx, sy = a.sum(), b.sum(); loo = (sx - a) / (sy - b); n = sel.sum()
            print('  %-12s offset %.2f slope %-6s %.3f +- %.3f' % (lab, off, nm, sx / sy, math.sqrt((n - 1) / n * np.sum((loo - loo.mean()) ** 2))))
