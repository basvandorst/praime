"""Puzzle 434: sigma sieve to N. Prints every x <= N/2 with sigma(sigma(x)) = sigma(x) + x
(needs sigma(x) <= N, true for every deficient x and so for every solution, see 434.md step 3),
every even x <= N with sigma(x) = 2x - 6, and the near misses |F| <= 100, F = sigma(sigma(x)) - sigma(x) - x,
per decade. Usage: p434_sieve.py 100000000 (about 15 s, 400 MB)."""
import math
import sys
import numpy as np

N = int(sys.argv[1])
S = math.isqrt(N)
s = np.zeros(N + 1, dtype=np.uint32)
for d in range(1, S + 1):
    s[d::d] += np.uint32(d)
for k in range(1, S + 1):
    hi = N // k
    if hi <= S:
        break
    s[k * (S + 1)::k] += np.arange(S + 1, hi + 1, dtype=np.uint32)

X = N // 2
x = np.arange(1, X + 1, dtype=np.int64)
y = s[1:X + 1].astype(np.int64)
ok = y <= N
x, y = x[ok], y[ok]
F = s[y].astype(np.int64) - y - x
print('solutions x <=', X, ':', x[F == 0].tolist())
even = np.arange(2, N + 1, 2, dtype=np.int64)
print('even x <=', N, 'with sigma(x) = 2x - 6:', even[2 * even - s[even].astype(np.int64) == 6].tolist())
composite = (y != x + 1) & (x > 1)
for lo, hi in ((10**4, 10**5), (10**5, 10**6), (10**6, 10**7), (10**7, X + 1)):
    m = composite & (np.abs(F) <= 100) & (x >= lo) & (x < hi)
    print('composite x in [%.0e, %.0e): %d near misses with |F| <= 100' % (lo, hi, m.sum()))

from sympy import divisor_sigma
extra = []
for m in range(1, math.isqrt(X) + 1):
    for t in (m * m, 2 * m * m):
        if t <= X:
            u = int(divisor_sigma(t))
            if u > N and int(divisor_sigma(u)) == u + t:
                extra.append(t)
print('squares and twice-squares x <= %d with sigma(x) > %d checked directly; solutions among them: %s' % (X, N, extra))
