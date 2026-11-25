"""stationary-phase prediction for a product's factor zeros under the product phase (post hoc, mechanism check):
sum_{rho zeros, gamma<=T} e^{2 i theta_K} ~ -c_rho * sum_{k <= T/2pi} rho(1/k)  (lines r = 1/k swept by the zeta chirp at t = 2 pi k),
with rho multiplicative, rho(p^-v) = p^{v/2} S_p(-v)/F_p(1/p)."""
import sys, os, math, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, cypari2
from ulib import *
pari = cypari2.Pari(); pari.allocatemem(300 * 10 ** 6)
for key, T in (('R23', 300), ('R31', 300)):
    spec = [c for c in CATALOGUE if c[0] == key][0][1]
    L = build(pari, spec); Q, k, mu = params(pari, L)
    K = int(T / (2 * math.pi))
    rho = {}
    for p in primes_upto(K):
        p = int(p); P = local_poly(pari, L, p, k); lam = series_inv(P, 80)
        F = sum(np.conj(lam[j]) * P[j] * p ** (-j) for j in range(len(P)))
        for v in range(1, 8):
            if p ** v > K: break
            s = sum(np.conj(lam[h - v]) * P[h] * p ** (-(2 * h - v) / 2) for h in range(v, len(P)))
            rho[(p, v)] = p ** (v / 2) * s / F
    tot = 0; part = []
    for kk in range(1, K + 1):
        f = pari.factor(kk); val = 1
        for i in range(len(f[0])):
            val *= rho[(int(f[0][i]), int(f[1][i]))]
        tot += val; part.append(tot)
    c = json.load(open('constants_%s.json' % key))['c']
    print('%s: K = %d, sum_{k<=K} rho(1/k) = %.3f%+.3fi, prediction -c*sum = %.2f ; partial sums at K/4, K/2: %.2f, %.2f' % (
        key, K, tot.real, tot.imag, -c * tot.real, part[K // 4 - 1].real, part[K // 2 - 1].real))
