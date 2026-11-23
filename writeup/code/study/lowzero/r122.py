"""Round 122 (preregister/round122.txt): Berry's per-prime rigidity at heights 900-1500 for 61 real characters.
usage: r122.py"""
import os, json, math, glob
import numpy as np
from sympy import primefactors
H = os.path.dirname(os.path.abspath(__file__)); ZD = os.path.join(H, '..', '..', '..', 'data', 'riemann', 'lzeros', 'r122')
def dS(p, L):
    return -4 / L ** 2 * sum(math.log(p) ** 2 * p ** (-k) * (1 - math.cos(2 * math.pi * k * math.log(p) / L)) / (2 * math.pi * k * math.log(p) / L) ** 2 for k in range(1, 40))
rows = []
for f in sorted(glob.glob(os.path.join(ZD, 'D*.npy'))):
    D = int(os.path.basename(f)[1:-4]); q = abs(D); g = np.load(f)
    u = g / (2 * math.pi) * np.log(q * g / (2 * math.pi * math.e))
    x = np.arange(u[0], u[-1] - 1, 0.01); cnt = np.searchsorted(u, x + 1) - np.searchsorted(u, x)
    L = float(np.mean(np.log(q * g / (2 * math.pi))))
    rows.append({'D': D, 'sigma2_1': float(np.var(cnt)), 'L': L, 'x': sum(dS(p, L) for p in primefactors(q)), 'nzeros': int(len(g))})
y = np.array([r['sigma2_1'] for r in rows]); X = np.array([[1, r['x']] for r in rows])
c = np.linalg.lstsq(X, y, rcond=None)[0]; res = y - X @ c; s2 = res @ res / (len(y) - 2); se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
d2_here = float(np.mean([dS(2, r['L']) for r in rows if r['D'] % 2 == 0]))
r117 = json.load(open(os.path.join(H, 'r117.json')))['pred']['2']
print('n %d, mean L %.2f, zeros per L %d; beta %.3f +- %.3f (t %.1f); predicted drop for 2|D here %.4f vs %.4f at 100-600' % (len(rows), np.mean([r['L'] for r in rows]), np.mean([r['nzeros'] for r in rows]), c[1], se[1], c[1] / se[1], d2_here, r117))
rd = [bool(c[1] > 3 * se[1]), bool(0.6 <= c[1] <= 1.5), bool(abs(d2_here) < abs(r117))]
print('readings', rd)
json.dump({'rows': rows, 'beta': float(c[1]), 'beta_se': float(se[1]), 'a': float(c[0]), 'pred_drop_2_here': d2_here, 'pred_drop_2_r117': r117, 'readings': rd}, open(os.path.join(H, 'r122.json'), 'w'), indent=0)
