"""v for Q = 210 at [1e9, 1e10) (n = 4.3e7 lags; the round-43 secondary that was not run), in chunks of k, with w(k)
by quadrature (Euler-Maclaurin) instead of an FFT and the prime tail by E1. Validated on [1e8, 1e9) against theory.F.
usage: v210_big.py L0 W"""
import sys, math, json
import numpy as np
from scipy.special import exp1
from scipy.interpolate import CubicSpline
import theory as T

Q = 210
L0 = int(float(sys.argv[1])); W = int(float(sys.argv[2]))
qf = [2, 3, 5, 7]
P = T.PR[(Q % T.PR) != 0]
base = float(np.sum(T.gen(P.astype(float)))) - 4 * exp1(math.log(T.PMAX))
c = math.prod(1 / (1 - 1 / p) ** 2 for p in qf) * math.exp(float(np.sum(np.log1p(-2.0 / P) - 2 * np.log1p(-1.0 / P))))
n = W // Q
rf = lambda x: c / (np.log(L0 + Q * x + Q / 2) * np.log(L0 + Q * x + Q / 2 + 2))
gx, gw = np.polynomial.legendre.leggauss(40)


def wsum(k):                                  # sum_{j=0}^{N} r(j) r(j+k), N = n-1-k, Euler-Maclaurin to first order
    N = n - 1 - k
    if N <= 0:
        return float(rf(0) * rf(k)) if N == 0 else 0.0
    edges = np.linspace(0, N, 21); tot = 0.0
    for a, b in zip(edges[:-1], edges[1:]):
        x = (b - a) / 2 * gx + (a + b) / 2
        tot += (b - a) / 2 * float(np.sum(gw * rf(x) * rf(x + k)))
    return tot + 0.5 * float(rf(0) * rf(k) + rf(N) * rf(N + k))


kg = np.unique(np.concatenate([np.linspace(1, n - 1, 3000).astype(np.int64), [1, 2, 3, n - 2, n - 1]]))
wg = np.array([wsum(int(k)) for k in kg])
ws = CubicSpline(kg, wg)
lim = int(math.isqrt(Q * n + 2)) + 1
Psmall = P[P <= lim]
F = 0.0; CH = 4 * 10 ** 6
for k0 in range(1, n, CH):
    k = np.arange(k0, min(k0 + CH, n), dtype=np.int64)
    d = Q * k
    L = np.full(k.size, base)
    for sh, nu in ((0, 2), (-2, 3), (2, 3)):
        cof = np.abs(d + sh) if sh else k.copy()          # p | Qk with p not | Q  <=>  p | k
        if sh:
            while True:
                m = cof % 2 == 0
                if not m.any(): break
                cof[m] //= 2
        for p in Psmall:
            p = int(p)
            j0 = ((-sh * pow(Q, -1, p)) % p - k0) % p
            idx = np.arange(j0, k.size, p)
            if idx.size == 0: continue
            L[idx] += math.log((1 - nu / p) / (1 - 2 / p) ** 2) - T.gen(p)
            while idx.size:
                cof[idx] //= p
                idx = idx[cof[idx] % p == 0]
        for q in qf:                                       # k may contain factors of Q: strip them
            while True:
                m = (cof % q == 0) & (cof > 1)
                if not m.any(): break
                cof[m] //= q
        big = cof > 1
        pb = cof[big].astype(float)
        L[big] += np.log1p(-nu / pb) - 2 * np.log1p(-2 / pb) - T.gen(pb)
    F += 2 * float(np.sum(ws(k) * (np.exp(L) - 1)))
xs = np.arange(n); mu = 0.0; s2 = 0.0
for a in range(0, n, CH):
    rr = rf(xs[a:a + CH].astype(float)); mu += float(rr.sum()); s2 += float(np.sum(rr * rr))
v = (mu - s2 + F) / mu
print(json.dumps(dict(L0=L0, W=W, n=n, v=v, mu=mu)))
