# Puzzle 997: checks every computed claim of results/997.md (about 10 seconds).
import numpy as np
from math import isqrt

# 1. Q1/Q2 (Volpatti): m - n + s - floor(m s / n) = 1 for all s < n < m  iff  (m - s + 1)^2 <= 4m + 1,  m < 600
for m in range(3, 600):
    for s in range(1, m):
        assert all(m - n + s - (m * s) // n == 1 for n in range(s + 1, m)) == ((m - s + 1)**2 <= 4 * m + 1)

# 2. primes below 2e8: the identity for (p_n, p_n+1, p_n+2) holds iff g_n g_n+1 < p_n+1, and (★) fails only for n = 3, 8
N = 200_000_001; sv = np.ones(N, bool); sv[:2] = False
for i in range(2, isqrt(N) + 1):
    if sv[i]: sv[i*i::i] = False
P = np.nonzero(sv)[0].astype(np.int64); del sv
g = np.diff(P); prod = g[:-1] * g[1:]; mid = P[1:-1]
fail = np.nonzero(prod >= mid)[0]
assert [(int(i) + 1, int(P[i])) for i in fail] == [(3, 5), (8, 19)]
for i in list(range(0, 200)) + list(range(0, len(P) - 2, 4999)):
    a, b, c = int(P[i]), int(P[i + 1]), int(P[i + 2])
    assert (c - b + a - (c * a) // b == 1) == ((b - a) * (c - b) < b)
assert int(g.max()) == 248                                   # largest gap below 2e8 (after 191912783)
big = np.nonzero(mid > 1000)[0]; r = prod[big] / mid[big]
assert int(P[big[np.argmax(r)]]) == 2161 and abs(float(r.max()) - 432 / 2179) < 1e-12   # 18 * 24 / 2179

# 3. above 2e8 and up to 4e18 every gap is at most 1476 (complete computation, Oliveira e Silva-Herzog-Pardi 2014;
#    OEIS A005250), so g_n g_n+1 <= 1476^2 = 2178576 < 2e8 < p_n+1
assert 1476**2 == 2178576 < 2 * 10**8
print("all claims verified")
