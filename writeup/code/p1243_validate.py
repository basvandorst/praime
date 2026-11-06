"""Puzzle 1243: runs of consecutive cube-root-smooth numbers (largest prime factor of k at most k^(1/3)).
Sieve of largest prime factors to N = 10^8: least start of a run of exactly n consecutive CRS numbers for n <= 5
(page: 8, 2400, 134848, 3678723, 84441507), the density of CRS numbers in [x/2, x] against Dickman's rho(3) = 0.0486,
and the number of maximal runs of each length against the independent model.
"""
import numpy as np
from math import log
N = 10**8
lpf = np.zeros(N + 1, dtype=np.int32)
for p in range(2, N + 1):
    if lpf[p] == 0:                       # p prime
        lpf[p::p] = p
k = np.arange(N + 1, dtype=np.int64)
crs = np.zeros(N + 1, dtype=bool)
small = lpf <= 465                                  # k <= 10^8 has cube root < 465, so larger lpf never qualifies
cube = np.where(small, lpf, 0).astype(np.int64)**3   # no overflow: 465^3 < 2^63
crs = small & (cube <= k)                            # lpf(k) <= k^(1/3)  <=>  lpf(k)^3 <= k
crs[0:2] = False
# maximal runs
d = np.diff(np.concatenate(([0], crs.astype(np.int8), [0])))
starts = np.flatnonzero(d == 1); ends = np.flatnonzero(d == -1)
lengths = ends - starts
least = {}
for s, L in zip(starts, lengths):
    if L not in least:
        least[int(L)] = int(s)
print("least start of a maximal run of exactly n CRS numbers:", dict(sorted(least.items())))
assert [least[n] for n in range(1, 6)] == [8, 2400, 134848, 3678723, 84441507]
for x in (10**6, 10**7, 10**8):
    dens = crs[x // 2:x + 1].mean()
    print(f"density of CRS numbers in [{x//2}, {x}]: {dens:.5f}   (Dickman rho(3) = 0.04861)")
counts = {int(L): int((lengths == L).sum()) for L in np.unique(lengths)}
print("number of maximal runs of each length below 10^8:", counts)
rho = 0.04861
print("independent model: runs of length n ~ N rho^n (1-rho)^2:", {n: round(N * rho**n * (1 - rho)**2) for n in range(1, 7)})
