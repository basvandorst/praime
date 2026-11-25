"""Partial sums S_H = sum_{h<=H} g(h), g multiplicative, g(p) = -2cos(2m theta_p)/(p-1) (p = 1 mod 4), g(p^2) = 0,
g = 0 on h with a factor 2 or q = 3 mod 4; compared with the Euler product (order of p) and c_exact*pi/4."""
import math, numpy as np
m = 4; H = 2_000_000
spf = np.zeros(H + 1, np.int64)
for q in range(2, H + 1):
    if spf[q] == 0:
        spf[q::q][spf[q::q] == 0] = q
# theta_p for p = 1 mod 4 <= H via a > b > 0
gp = np.zeros(H + 1)
R = math.isqrt(H) + 1
for a in range(1, R + 1):
    b = np.arange(1, a)
    p = a * a + b * b
    ok = (p <= H)
    p = p[ok]; bb = b[ok]
    isp = spf[p] == p
    gp[p[isp]] = -2 * np.cos(2 * m * np.arctan2(bb[isp], a)) / (p[isp] - 1)
g = np.zeros(H + 1); g[1] = 1.0
for h in range(2, H + 1):
    p = spf[h]; k = h // p
    g[h] = 0.0 if k % p == 0 else g[k] * gp[p]
S = np.cumsum(g)
c_exact = 2.026446079853554
target = c_exact * math.pi / 4
print("target (pi/4) c_exact =", target)
prod = 1.0
pp = np.nonzero(gp)[0]
lp = np.cumsum(np.log1p(gp[pp]))
for Hx in (10**2, 10**3, 10**4, 10**5, 10**6, 2 * 10**6):
    hh = np.arange(1, Hx + 1)
    ces = float(np.sum(g[1:Hx + 1] * (1 - hh / Hx)))
    ces2 = float(np.sum(g[1:Hx + 1] * (1 - hh / Hx) ** 2))
    ip = np.searchsorted(pp, Hx, side='right')
    print(f"H={Hx:>8d}  S_H={S[Hx]:.6f}  Riesz1={ces:.6f}  Riesz2={ces2:.6f}  prod_p<=H={math.exp(lp[ip-1]):.6f}   S_H-target={S[Hx]-target:+.2e}")
