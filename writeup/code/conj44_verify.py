# Conjecture 44: for consecutive primes p < q, z = sqrt((p^2 + q^2)/2 - 1) is an integer iff q - p = 2.  (a few seconds)
import numpy as np
from math import isqrt
N = 200_000_001; sv = np.ones(N, bool); sv[:2] = False
for i in range(2, isqrt(N) + 1):
    if sv[i]: sv[i*i::i] = False
P = np.nonzero(sv)[0].astype(np.int64); del sv
p, q = P[1:-1], P[2:]                                             # odd consecutive primes below 2e8
g = q - p
# 1. Oakes' reduction: with q = p + 2t, z^2 = (p + t)^2 + t^2 - 1, strictly between (p+t)^2 and (p+t+1)^2 when 1 < t < sqrt(2p+3) + 1
t = g // 2
assert np.all((t == 1) | ((t - 1)**2 < 2 * p + 3))                # no gap below 2e8 reaches the bound
# 2. direct check of the statement itself on every pair below 2e8 (exact integer square test)
z2 = (p * p + q * q) // 2 - 1                                      # p, q odd: p^2 + q^2 = 2 (mod 4), so the halving is exact
zr = np.sqrt(z2.astype(np.float64)).astype(np.int64)
issq = np.zeros(len(z2), bool)                                    # exact square test (float root +- 1)
for d in (-1, 0, 1):
    zz = zr + d; issq |= zz * zz == z2
assert np.array_equal(issq, g == 2)
# 3. above 2e8 and up to 4e18 every gap is <= 1476, so t <= 738 and (t - 1)^2 <= 737^2 < 2 p + 3
assert 737**2 < 2 * 2 * 10**8 + 3
print("all claims verified")
