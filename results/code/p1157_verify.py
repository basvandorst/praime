# Puzzle 1157: for consecutive primes p < q,
#   (1) floor((p^2+q^2)/(2p)) = q  <=> g^2 < 2p,   (2) floor((p^2+q^2)/(2q)) = p <=> g^2 < 2q,
#   (3) floor((p^2+q^2)/(p+q)) = (p+q)/2 <=> g^2 < 2(p+q)   (p > 2),   g = q - p.   (a few seconds)
import numpy as np
from math import isqrt
N = 200_000_001; sv = np.ones(N, bool); sv[:2] = False
for i in range(2, isqrt(N) + 1):
    if sv[i]: sv[i*i::i] = False
P = np.nonzero(sv)[0].astype(np.int64); del sv
p, q = P[:-1], P[1:]; g = q - p; s = p * p + q * q
f1 = s // (2 * p) == q; f2 = s // (2 * q) == p
odd = p > 2; f3 = (s[odd] // (p[odd] + q[odd])) == (p[odd] + q[odd]) // 2
assert np.array_equal(f1, g * g < 2 * p) and np.array_equal(f2, g * g < 2 * q) and np.array_equal(f3, (g * g < 2 * (p + q))[odd])
assert [int(x) for x in p[~f1]] == [7] and not (~f2).any() and not (~f3).any()     # the only failure below 2e8: formula 1 at p = 7
# above 2e8 (to 4e18) every gap is <= 1476 and 1476^2 = 2178576 < 2p
assert 1476**2 < 2 * 2 * 10**8
print("all claims verified")
