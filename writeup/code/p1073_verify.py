"""Puzzle 1073: p*q = 1 + k^2 has no solution in consecutive primes p < q with p <= 4e18.  (about 3 s)

1. Reduction (algebra, checked numerically): a solution needs q - p >= 2 + 2*sqrt(2p - 1).
2. Direct check: no consecutive pair below 1e7 gives a square p*q - 1, and every gap there is below the bound.
3. Above 1e7 the bound exceeds 1476, the largest gap below 4e18 (complete computation, literature).
4. The bound is sharp: e with 2e^2 - 2e + 1 and 2e^2 + 2e + 1 both prime gives p*q = 1 + (2e^2)^2, q - p = 4e.
"""
import numpy as np
from math import isqrt
from sympy import isprime

N = 10**7
sv = np.ones(N + 200, bool); sv[:2] = False
for i in range(2, isqrt(N + 200) + 1):
    if sv[i]: sv[i*i::i] = False
P = np.nonzero(sv)[0].astype(np.int64)
p, q = P[:-1], P[1:]
keep = p < N; p, q = p[keep], q[keep]
v = p * q - 1                                        # < 1.0001e14, exact in int64
r = np.sqrt(v.astype(np.float64)).astype(np.int64)
sq = np.zeros(len(v), bool)
for d in (-1, 0, 1):
    sq |= (r + d) * (r + d) == v
assert not sq.any()                                  # no consecutive pair below 1e7
g = q - p
assert np.all((g - 2)**2 < 4 * (2 * p - 1))          # every gap below the bound 2 + 2 sqrt(2p - 1)
# the reduction itself, on all pairs p < q < 2000 (not only consecutive): a solution has q - p >= 2 + 2 sqrt(2p - 1)
for a in range(3, 2000, 2):
    for b in range(a + 2, 2000, 2):
        k2 = a * b - 1
        if isqrt(k2)**2 == k2: assert (b - a - 2)**2 >= 4 * (2 * a - 1), (a, b)
# above 1e7 the bound is > 2 + 2*sqrt(2e7 - 1) = 8946 > 1476, the largest gap below 4e18
assert (1476 - 2)**2 < 4 * (2 * N - 1)
# sharpness: Sophie Germain identity 4e^4 + 1 = (2e^2 - 2e + 1)(2e^2 + 2e + 1)
fam = []
for e in range(1, 200):
    a, b = 2*e*e - 2*e + 1, 2*e*e + 2*e + 1
    assert a * b == 1 + (2*e*e)**2 and (b - a - 2)**2 == 4 * (2 * a - 1)
    if isprime(a) and isprime(b): fam.append((e, a, b))
print("Germain family with both prime (e, p, q):", fam[:8])
print("all claims verified")
