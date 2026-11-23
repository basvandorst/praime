"""Round 128 (preregister/round128.txt): zeta vs real characters of prime conductor at equal density L in [7.6, 8.6]:
Sigma2(1) difference against Berry's dS(q, L). usage: r128.py"""
import os, json, math, glob
import numpy as np
H = os.path.dirname(os.path.abspath(__file__)); RD = os.path.join(H, '..', '..', '..', 'data', 'riemann')
def dS(p, L): return -(2 / math.pi ** 2) * sum(math.sin(math.pi * k * math.log(p) / L) ** 2 / (k * k * p ** k) for k in range(1, 60))
def stats(u):
    x = np.arange(u[0], u[-1] - 1, 0.01); cnt = np.searchsorted(u, x + 1) - np.searchsorted(u, x)
    s2 = float(np.var(cnt)); sp = np.diff(u)
    nb = len(u) // 200; bv = []
    for b in range(nb):
        ub = u[b * 200:(b + 1) * 200]; xb = np.arange(ub[0], ub[-1] - 1, 0.01); cb = np.searchsorted(ub, xb + 1) - np.searchsorted(ub, xb); bv.append(np.var(cb))
    return s2, float(np.std(bv) / math.sqrt(nb)), float(np.var(sp))
g = np.array([float(l) for l in open(os.path.join(RD, 'zeros_40000_flint.txt')) if l.strip()])
L = np.log(g / (2 * math.pi)); m = (L >= 7.6) & (L <= 8.6)
uz = (g / (2 * math.pi) * np.log(g / (2 * math.pi * math.e)) + 7 / 8)[m]
z2, zse, zsp = stats(uz)
rows = []
for f in sorted(glob.glob(os.path.join(RD, 'lzeros', 'r128', 'D*.npy'))):
    D = int(os.path.basename(f)[1:-4]); q = abs(D); t = np.load(f)
    u = t / (2 * math.pi) * np.log(q * t / (2 * math.pi * math.e))
    s2, se, sp = stats(u); Lm = float(np.mean(np.log(q * t / (2 * math.pi))))
    rows.append({'D': D, 'sigma2_1': s2, 'se': se, 'spacing_var': sp, 'L': Lm, 'pred': dS(q, Lm), 'n': int(len(t))})
diff = np.mean([r['sigma2_1'] for r in rows]) - z2; pred = np.mean([r['pred'] for r in rows])
se = math.sqrt(np.mean([r['se'] ** 2 for r in rows]) / len(rows) + zse ** 2)
r1 = abs(diff - pred) <= 2 * se
print('zeta: Sigma2(1) %.4f +- %.4f (n %d, L mean %.2f), spacing var %.4f' % (z2, zse, int(m.sum()), float(L[m].mean()), zsp))
for r in rows: print('D %4d n %4d L %.2f Sigma2 %.4f +- %.4f pred diff %+.4f spacing %.4f' % (r['D'], r['n'], r['L'], r['sigma2_1'], r['se'], r['pred'], r['spacing_var']))
print('mean difference chi - zeta %+.4f +- %.4f, predicted %+.4f; reading (1): %s; spacing-variance difference %+.4f' % (diff, se, pred, r1, np.mean([r['spacing_var'] for r in rows]) - zsp))
json.dump({'zeta': {'sigma2_1': z2, 'se': zse, 'spacing_var': zsp}, 'rows': rows, 'mean_diff': float(diff), 'se': se, 'pred': float(pred), 'reading1': bool(r1)}, open(os.path.join(H, 'r128.json'), 'w'), indent=1)
