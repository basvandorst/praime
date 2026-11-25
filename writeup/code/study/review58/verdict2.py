import json, math, sys, glob
import numpy as np
sys.path.insert(0, '.')
from r1_offsets import theta, jk_mean, load
old = [r for r in load() if r['cls'] == 'g2']
def res(g, Q, mu, eps, c):
    gg = g[(g > 0) & (g <= 150)]; return np.exp(2j * theta(gg, Q, mu)).sum().real + eps * c * 150 / (2 * math.pi)
xo = np.array([res(r['g'], r['Q'], r['mu'], r['eps'], r['c']) for r in old]); No = np.array([r['Q'] for r in old]); eo = np.array([r['eps'] for r in old])
rep = []
for fn in sorted(glob.glob('rep/*.npz')):
    z = np.load(fn); i = json.loads(str(z['info'])); rep.append((res(z['g'], i['Q'], i['mu'], i['eps'][0], i['c']), i['Q'], i['eps'][0]))
xr = np.array([a for a, _, _ in rep]); Nr = np.array([b for _, b, _ in rep]); er = np.array([e for _, _, e in rep])
for lab, sel in (('old N < 1055', No < 1055), ('old N >= 1055', No >= 1055), ('old N < 2000', No < 2000), ('old N >= 2000', No >= 2000)):
    m, s = jk_mean(xo[sel]); print('%-15s n=%2d  %+.2f +- %.2f   eps+1 share %.2f' % (lab, sel.sum(), m, s, (eo[sel] > 0).mean()))
for lab, sel in (('rep N < 2000', Nr < 2000), ('rep N >= 2000', Nr >= 2000)):
    m, s = jk_mean(xr[sel]); print('%-15s n=%2d  %+.2f +- %.2f   eps+1 share %.2f' % (lab, sel.sum(), m, s, (er[sel] > 0).mean()))
xa = np.concatenate([xo, xr]); Na = np.concatenate([No, Nr]); src = np.r_[np.ones(len(xo)), np.zeros(len(xr))]
sel = Na >= 1055; m, s = jk_mean(xa[sel & (src == 1)]); m2, s2 = jk_mean(xa[sel & (src == 0)])
print('same range N >= 1055: old %+.2f +- %.2f (n=%d), rep %+.2f +- %.2f (n=%d), diff %.2f +- %.2f' % (m, s, (sel & (src == 1)).sum(), m2, s2, (sel & (src == 0)).sum(), m - m2, math.hypot(s, s2)))
# regression of residual on log N, all 100, with a source dummy
X = np.c_[np.ones(len(xa)), np.log(Na) - np.log(Na).mean(), src]; b, *_ = np.linalg.lstsq(X, xa, rcond=None)
r = xa - X @ b; cov = np.linalg.inv(X.T @ X) * (r @ r) / (len(xa) - 3); print('fit res = a + b log N + d*old: a %.2f b %.2f +- %.2f d %.2f +- %.2f' % (b[0], b[1], math.sqrt(cov[1, 1]), b[2], math.sqrt(cov[2, 2])))
