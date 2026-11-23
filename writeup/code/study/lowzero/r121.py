"""Round 121 (preregister/round121.txt): Sigma2(1) of complex characters against x = sum_{p|N} dS(p, L); fit
Sigma2(1) = a + beta x (Berry's diagonal term: beta = 1). usage: r121.py"""
import os, json, math
import numpy as np
from sympy import primefactors
H = os.path.dirname(os.path.abspath(__file__)); ZD = os.path.join(H, '..', '..', '..', 'data', 'riemann', 'lzeros')
Z = {}
for i in range(4):
    z = np.load(os.path.join(ZD, 'r121_part%d.npz' % i)); Z.update({k: z[k] for k in z.files})
def dS(p, L):
    return -4 / L ** 2 * sum(math.log(p) ** 2 * p ** (-k) * (1 - math.cos(2 * math.pi * k * math.log(p) / L)) / (2 * math.pi * k * math.log(p) / L) ** 2 for k in range(1, 40))
rows = []
for key, z in Z.items():
    N = int(key.split('_')[0][1:])
    vs, ws, Ls = [], [], []
    for side in (1, -1):
        g = np.sort(side * z[side * z > 0])
        u = g / (2 * math.pi) * np.log(N * g / (2 * math.pi * math.e))
        x = np.arange(u[0], u[-1] - 1, 0.01)
        cnt = np.searchsorted(u, x + 1) - np.searchsorted(u, x)
        vs.append(np.var(cnt)); ws.append(len(x)); Ls.append(np.mean(np.log(N * g / (2 * math.pi))))
    L = float(np.average(Ls, weights=ws)); s2 = float(np.average(vs, weights=ws))
    xpred = sum(dS(p, L) for p in primefactors(N))
    rows.append({'key': key, 'N': N, 'sigma2_1': s2, 'L': L, 'x': xpred, 'primes': [int(p) for p in primefactors(N)]})
y = np.array([r['sigma2_1'] for r in rows]); X = np.array([[1, r['x']] for r in rows])
c = np.linalg.lstsq(X, y, rcond=None)[0]; res = y - X @ c; s2 = res @ res / (len(y) - 2); se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
beta, bse = float(c[1]), float(se[1])
groups = {}
for r in rows:
    g = 'prime' if len(r['primes']) == 1 else ('%d|N' % min(r['primes']))
    groups.setdefault(g, []).append(r['sigma2_1'])
print({g: (len(v), round(float(np.mean(v)), 4)) for g, v in groups.items()})
print('n %d, beta %.3f +- %.3f (t = %.1f), a %.4f; readings: %s %s' % (len(rows), beta, bse, beta / bse, c[0], beta > 4 * bse, 0.7 <= beta <= 1.4))
json.dump({'rows': rows, 'beta': beta, 'beta_se': bse, 'a': float(c[0]), 'groups': {g: [len(v), float(np.mean(v))] for g, v in groups.items()}, 'readings': [bool(beta > 4 * bse), bool(0.7 <= beta <= 1.4)]}, open(os.path.join(H, 'r121.json'), 'w'), indent=0)
