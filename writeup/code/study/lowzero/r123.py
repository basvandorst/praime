"""Round 123 (preregister/round123.txt): Sigma2(l) for l = 0.5 .. 8 on the 181 real L-functions, per-prime fits against
Berry's oscillating prediction -(2/pi^2) sum_k sin^2(pi l k ln p/L)/(k^2 p^k). usage: r123.py"""
import os, json, math
import numpy as np
H = os.path.dirname(os.path.abspath(__file__)); ZD = os.path.join(H, '..', '..', '..', 'data', 'riemann', 'lzeros')
Z = {}
for f in ('r117_part0.npz', 'r117_part1.npz', 'r120_part0.npz', 'r120_part1.npz'):
    z = np.load(os.path.join(ZD, f)); Z.update({k: z[k] for k in z.files})
PR = (2, 3, 5, 7, 11, 13); LS = [0.5 * i for i in range(1, 17)]
def dS(p, L, l): return -(2 / math.pi ** 2) * sum(math.sin(math.pi * l * k * math.log(p) / L) ** 2 / (k * k * p ** k) for k in range(1, 60))
rows = []
for key, g in Z.items():
    D = int(key[1:]); q = abs(D); m = (g >= 100) & (g <= 600)
    u = (g / (2 * math.pi) * np.log(q * g / (2 * math.pi * math.e)))[m]
    s2 = {}
    for l in LS:
        x = np.arange(u[0], u[-1] - l, 0.01); cnt = np.searchsorted(u, x + l) - np.searchsorted(u, x); s2[l] = float(np.var(cnt))
    rows.append({'D': D, 'L': float(np.mean(np.log(q * g[m] / (2 * math.pi)))), 's2': s2})
X = np.array([[1] + [int(r['D'] % p == 0) for p in PR] for r in rows], float); XtXi = np.linalg.inv(X.T @ X)
fit, pred, se = {}, {}, {}
for l in LS:
    y = np.array([r['s2'][l] for r in rows]); c = XtXi @ X.T @ y; res = y - X @ c; s = math.sqrt(res @ res / (len(y) - X.shape[1]))
    for i, p in enumerate(PR):
        fit[(p, l)] = float(c[i + 1]); se[(p, l)] = float(s * math.sqrt(XtXi[i + 1, i + 1]))
        pred[(p, l)] = float(np.mean([dS(p, r['L'], l) for r in rows if r['D'] % p == 0]))
f = np.array([fit[k] for k in fit]); pv = np.array([pred[k] for k in fit]); slope = float(f @ pv / (pv @ pv))
print('slope through origin over 96 points: %.3f' % slope)
for p in PR: print('p=%2d' % p, ' '.join('%+.3f/%+.3f' % (fit[(p, l)], pred[(p, l)]) for l in (1, 2, 3.5, 5.5, 6.5, 7.5)))
r1 = 0.7 <= slope <= 1.3; r2 = abs(fit[(2, 6.5)]) > 5 * abs(fit[(2, 1.0)]); r3 = abs(fit[(13, 3.5)]) < abs(fit[(13, 2.0)]) and abs(fit[(13, 7.5)]) > abs(fit[(13, 3.5)])
print('readings', r1, r2, r3, ' p13: l=2 %+.4f l=3.5 %+.4f l=7.5 %+.4f (se ~%.4f)' % (fit[(13, 2.0)], fit[(13, 3.5)], fit[(13, 7.5)], se[(13, 3.5)]))
json.dump({'slope': slope, 'fit': {'%d_%s' % k: v for k, v in fit.items()}, 'pred': {'%d_%s' % k: v for k, v in pred.items()}, 'se': {'%d_%s' % k: v for k, v in se.items()}, 'readings': [bool(r1), bool(r2), bool(r3)]}, open(os.path.join(H, 'r123.json'), 'w'), indent=0)
