# Puzzle 211 (attempt): AS(p_n) = p_n - p_(n-1) + ... = a(n), a(n) = p_n - a(n-1).  Numerical behaviour of a(n)/p_n
# and of the bound used in writeup/211.md: the share of gaps of one index parity, for primes up to 2e8.
import numpy as np
from math import isqrt
N = 200_000_001; sv = np.ones(N, bool); sv[:2] = False
for i in range(2, isqrt(N) + 1):
    if sv[i]: sv[i*i::i] = False
P = np.nonzero(sv)[0].astype(np.int64)
a = np.empty(len(P), np.int64); s = 0
for i, p in enumerate(P.tolist()):
    s = p - s; a[i] = s
b = a / P
assert 0.4999 < b[10**7 - 1] < 0.5001 and np.all(np.abs(b[10**6:] - 0.5) < 0.001)
assert np.max(np.abs(2 * a[1000:] - P[1000:]) / np.sqrt(P[1000:])) < 8
# sieve-bound ingredient: fraction of primes p <= X whose next prime is within (1/16) ln X
g = np.diff(P); X = float(P[-1]); frac = np.mean(g <= np.log(X) / 16)
assert frac < 4 / 16                       # the theorem allows up to 4 * delta = 1/4 asymptotically; real value is far smaller
print("b(n) in", b[10**6:].min(), b[10**6:].max(), "for 1e6 <= n <=", len(P), "; short-gap fraction", round(frac, 4))
print("all claims verified")
