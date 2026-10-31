# Puzzle 720: (m+-1)^2 * 10^k + m^2 = Z^2 with k = number of digits of m^2 (m = right block, m +- 1 = left block).
from math import isqrt
sols = []
for k in range(1, 16):
    lo = isqrt(10**(k - 1) - 1) + 1 if k > 1 else 1      # smallest m with m^2 >= 10^(k-1)
    hi = isqrt(10**k - 1)                                 # largest m with m^2 < 10^k
    P = 10**k
    for m in range(max(lo, 1), hi + 1):
        r = m * m
        for s in (m - 1, m + 1):
            if s < 1: continue
            N = s * s * P + r
            z = isqrt(N)
            if z * z == N: sols.append((k, s, m, N))
    print(f"k={k:2d}: m in [{lo}, {hi}] checked, solutions so far: {sols}", flush=True)
assert all(min(s, m) <= 3 for k, s, m, N in sols), sols
print("no solution with x, y > 3 for k <= 15 (every m below 10^7.5)")
