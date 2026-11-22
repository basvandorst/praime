"""Round 120 (preregister/round120.txt): Sigma2(1) on a dense grid for all fundamental discriminants 100 <= |D| <= 400,
fitted on indicators of p | D, p = 2, 3, 5, 7, 11, 13, against Berry's per-prime prediction. usage: r120.py"""
import os, json, math
import numpy as np
H = os.path.dirname(os.path.abspath(__file__)); ZD = os.path.join(H, '..', '..', '..', 'data', 'riemann', 'lzeros')
Z = {}
for f in ('r117_part0.npz', 'r117_part1.npz', 'r120_part0.npz', 'r120_part1.npz'):
    z = np.load(os.path.join(ZD, f)); Z.update({k: z[k] for k in z.files})
PR = (2, 3, 5, 7, 11, 13)
def dS(p, L):
    return -4 / L ** 2 * sum(math.log(p) ** 2 * p ** (-k) * (1 - math.cos(2 * math.pi * k * math.log(p) / L)) / (2 * math.pi * k * math.log(p) / L) ** 2 for k in range(1, 40))
rows = []
for k, g in Z.items():
    D = int(k[1:]); q = abs(D)
    u = g / (2 * math.pi) * np.log(q * g / (2 * math.pi * math.e)); m = (g >= 100) & (g <= 600); uu = u[m]
    x = np.arange(uu[0], uu[-1] - 1, 0.01)
    cnt = np.searchsorted(uu, x + 1) - np.searchsorted(uu, x)
    rows.append({'D': D, 'sigma2_1': float(np.var(cnt)), 'L': float(np.mean(np.log(q * g[m] / (2 * math.pi)))), **{'div%d' % p: int(D % p == 0) for p in PR}})
y = np.array([r['sigma2_1'] for r in rows]); X = np.array([[1] + [r['div%d' % p] for p in PR] for r in rows], float)
c = np.linalg.lstsq(X, y, rcond=None)[0]; res = y - X @ c; s2 = res @ res / (len(y) - X.shape[1]); cov = s2 * np.linalg.inv(X.T @ X); se = np.sqrt(np.diag(cov))
pred = {p: float(np.mean([dS(p, r['L']) for r in rows if r['div%d' % p]])) for p in PR}
fit = {p: (float(c[i + 1]), float(se[i + 1])) for i, p in enumerate(PR)}
for p in PR: print('p=%2d n=%3d fitted %+.4f +- %.4f  predicted %+.4f' % (p, sum(r['div%d' % p] for r in rows), fit[p][0], fit[p][1], pred[p]))
f = np.array([fit[p][0] for p in PR]); pv = np.array([pred[p] for p in PR]); slope = float(f @ pv / (pv @ pv))
idx = [PR.index(p) + 1 for p in (7, 11, 13)]; w = np.zeros(len(c)); w[idx] = 1 / 3
pool = float(w @ c); pool_se = float(math.sqrt(w @ cov @ w))
r1 = int((f < 0).sum()) >= 5; r2 = 0.7 <= slope <= 1.6; r3 = pool < -2 * pool_se
print('n %d; slope through origin %.2f; pooled 7/11/13: %+.4f +- %.4f; readings %s %s %s' % (len(rows), slope, pool, pool_se, r1, r2, r3))
json.dump({'rows': rows, 'fit': fit, 'pred': pred, 'slope': slope, 'pooled_7_11_13': [pool, pool_se], 'readings': [r1, r2, r3]}, open(os.path.join(H, 'r120.json'), 'w'), indent=0)
