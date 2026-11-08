"""Conjecture 60 (Noe): for n >= 1 and k > K = log 127/log 16 there is a prime strictly between n^k and (n+1)^k.
Proof for every n, k with n^k <= 4e18.  (about 2 s)

A failure at (n, k) means consecutive primes p < p' with p <= n^k and (n+1)^k <= p', i.e.
    ln p / ln n <= k <= ln p' / ln(n+1).
Then p' - p >= (n+1)^k - n^k >= k n^(k-1) >= k p^(1-1/k) >= K p^(1-1/K)   (k p^(1-1/k) increases with k).
Below 4e18 every gap is <= 1476, so p < (1476/K)^(K/(K-1)) < 7.0e6.  For every pair of consecutive primes
below 7.0e6 with p' - p >= K p^(1-1/K) we list all n with a nonempty interval of failing k, and check that
none reaches beyond K; the largest failing k is exactly K (n = 15, gap 113 -> 127, 16^K = 127).
"""
import math
import numpy as np
from mpmath import mp, mpf, log as mlog
mp.dps = 50
K = math.log(127) / math.log(16)
PB = (1476 / K) ** (K / (K - 1))
print(f"K = {K:.12f}; failures need a gap >= K p^(1-1/K); with gaps <= 1476 this needs p < {PB:.4e}")
N = int(PB) + 2000
sv = np.ones(N, bool); sv[:2] = False
for i in range(2, int(N**0.5) + 1):
    if sv[i]: sv[i*i::i] = False
P = np.nonzero(sv)[0]
cand = [(int(a), int(b)) for a, b in zip(P[:-1], P[1:]) if a < PB and b - a >= K * a ** (1 - 1 / K) * (1 - 1e-9)]
print(f"{len(cand)} prime pairs below {PB:.3e} have a gap that large; the largest p among them is {max(a for a, b in cand)}")
worst, nfail = None, 0
for p, q in cand:
    for n in range(2, int(q ** (1 / K)) + 2):
        lo, hi = math.log(p) / math.log(n), math.log(q) / math.log(n + 1)
        if lo <= hi:                                   # failing k in [lo, hi]
            nfail += 1
            if worst is None or hi > worst[0]: worst = (hi, lo, n, p, q)
hi, lo, n, p, q = worst
print(f"{nfail} failing intervals with k >= 1 from these pairs; the highest ends at k = {hi:.12f} (n = {n}, gap {p} -> {q})")
assert (n, p, q) == (15, 113, 127)
assert mlog(q) / mlog(n + 1) == mlog(127) / mlog(16)   # exactly K: the open interval (15^K, 127) has no prime
# every other failing interval ends strictly below K (checked with 50-digit arithmetic for near ties)
for p, q in cand:
    for m in range(2, int(q ** (1 / K)) + 2):
        a, b = mlog(p) / mlog(m), mlog(q) / mlog(m + 1)
        if a <= b and (m, p, q) != (15, 113, 127):
            assert b < mlog(127) / mlog(16), (m, p, q)
print("no failure with k > K and n^k <= 4e18; at k = K exactly the statement fails at n = 15")
print("all claims verified")
