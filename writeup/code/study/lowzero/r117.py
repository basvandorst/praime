"""Round 117 (preregister/round117.txt): number variance Sigma2(1) of the unfolded zeros (100 <= gamma <= 600) of the 121
L-functions of rounds 115-116, fitted on indicators of 2|D, 3|D, 5|D, against Berry's diagonal prediction
dSigma2(1) = -(4/L^2) sum_k (ln p)^2 p^-k (1 - cos a_k)/a_k^2, a_k = 2 pi k ln p / L. usage: r117.py"""
import os, json, math
import numpy as np
H = os.path.dirname(os.path.abspath(__file__)); ZD = os.path.join(H, '..', '..', '..', 'data', 'riemann', 'lzeros')
Z = {}
for f in ('r117_part0.npz', 'r117_part1.npz'):
    z = np.load(os.path.join(ZD, f)); Z.update({k: z[k] for k in z.files})
rng = np.random.default_rng(117); rows = []
for k, g in Z.items():
    D = int(k[1:]); q = abs(D)
    u = g / (2 * math.pi) * np.log(q * g / (2 * math.pi * math.e))
    m = (g >= 100) & (g <= 600); uu = u[m]
    x = rng.uniform(uu[0], uu[-1] - 1, 40000)
    cnt = np.searchsorted(uu, x + 1) - np.searchsorted(uu, x)
    Lmean = float(np.mean(np.log(q * g[m] / (2 * math.pi))))
    rows.append({'D': D, 'sigma2_1': float(np.var(cnt)), 'L': Lmean, 'div2': int(D % 2 == 0), 'div3': int(D % 3 == 0), 'div5': int(D % 5 == 0)})
def dS(p, L):
    tot = 0.0
    for kk in range(1, 40):
        a = 2 * math.pi * kk * math.log(p) / L
        tot += math.log(p) ** 2 * p ** (-kk) * (1 - math.cos(a)) / a ** 2
    return -4 / L ** 2 * tot
Ls = [r['L'] for r in rows]
pred = {p: float(np.mean([dS(p, L) for L in Ls])) for p in (2, 3, 5)}
y = np.array([r['sigma2_1'] for r in rows]); X = np.array([[1, r['div2'], r['div3'], r['div5']] for r in rows], float)
c, res, rk, sv = np.linalg.lstsq(X, y, rcond=None)
resid = y - X @ c; s2 = resid @ resid / (len(y) - 4); cov = s2 * np.linalg.inv(X.T @ X); se = np.sqrt(np.diag(cov))
fit = {p: (float(c[i]), float(se[i])) for i, p in ((1, 2), (2, 3), (3, 5))}
for p in (2, 3, 5): print('p=%d fitted %+.4f +- %.4f  predicted %+.4f  ratio %.2f' % (p, fit[p][0], fit[p][1], pred[p], fit[p][0] / pred[p]))
print('mean L %.2f, constant %.4f, n %d' % (np.mean(Ls), c[0], len(rows)))
r1 = all(fit[p][0] < 0 for p in (2, 3, 5)); r2 = all(0.5 <= fit[p][0] / pred[p] <= 2 for p in (2, 3, 5))
r3 = sorted((2, 3, 5), key=lambda p: fit[p][0]) == sorted((2, 3, 5), key=lambda p: pred[p])
print('readings:', r1, r2, r3)
json.dump({'rows': rows, 'fit': fit, 'pred': pred, 'const': float(c[0]), 'readings': [r1, r2, r3]}, open(os.path.join(H, 'r117.json'), 'w'), indent=0)
