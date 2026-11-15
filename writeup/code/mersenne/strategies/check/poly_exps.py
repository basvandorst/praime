"""Exponents on the corner rays of the polygon spirals (P = 3..40, rays c_j(m) = P m(m-1)/2 + j m + 1, m >= 1, <= 500,000).
Expected count under the heuristic w(p) = e^gamma log2(a p)/p (a = 6 for p = 1 mod 4, 2 for p = 3 mod 4; w = 1 for p <= 7), Poisson tail;
family-wise: the largest exponent count over all 2,698 rays vs 2,000 simulated exponent sets drawn with weights w(p)."""
import numpy as np, math, sympy
N = 500000
E = [p for p in [2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243,110503,132049,216091] if p <= N]
s = np.ones(N + 1, bool); s[:2] = False
for i in range(2, int(N ** .5) + 1):
    if s[i]: s[i * i::i] = False
pr = np.flatnonzero(s)
w = np.array([1.0 if p <= 7 else min(1.0, math.exp(0.5772156649) * math.log2((6 if p % 4 == 1 else 2) * p) / p) for p in pr])
def pois_tail(o, e):
    return 1 - sum(math.exp(-e) * e ** i / math.factorial(i) for i in range(o))
rays = []
for P in range(3, 41):
    for j in range(1, P + 1):
        m = np.arange(1, 2000); f = P * m * (m - 1) // 2 + j * m + 1; f = f[f <= N]
        rays.append((P, j, f))
Eset = set(E); widx = {p: i for i, p in enumerate(pr)}
best = []
for P, j, f in rays:
    o = sum(1 for v in f if v in Eset); e = sum(w[widx[v]] for v in f if s[v])
    best.append((o, e, P, j, [int(v) for v in f if v in Eset]))
b6 = [b for b in best if b[2] == 6 and b[3] == 6][0]
print(f'P=6, ray 3m^2+3m+1: exponents {b6[4]}; expected {b6[1]:.2f}; Poisson P(>= {b6[0]}) = {pois_tail(b6[0], b6[1]):.3f}')
top = sorted(best, key=lambda b: pois_tail(b[0], b[1]))[:5]
for o, e, P, j, ex in top: print(f'  P={P} j={j}: {o} exponents {ex}, expected {e:.2f}, tail {pois_tail(o, e):.4f}')
rng = np.random.default_rng(7); pmin_real = min(pois_tail(b[0], b[1]) for b in best)
ray_sets = [set(int(v) for v in f) for _, _, f in rays]; ray_e = [b[1] for b in best]
sims = []
for _ in range(2000):
    lab = set(rng.choice(pr, size=len(E), replace=False, p=w / w.sum()).tolist())
    sims.append(min(pois_tail(len(lab & rs), e) for rs, e in zip(ray_sets, ray_e)))
sims = np.array(sims)
print(f'smallest tail over all {len(rays)} rays: {pmin_real:.4f}; family-wise p (2,000 simulated exponent sets) = {(1 + (sims <= pmin_real).sum()) / 2001:.3f}')
