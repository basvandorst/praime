"""Puzzle 309: f(x) = sigma(phi(x)) - phi(sigma(x)) for every x <= N (sieves: phi to 6N, sigma to N).
Prints the solutions, near misses |f(x) - x| <= 100 per decade, the counts of F = f(x) - x near 0 for
even x, and for each class x mod 24 how often f(x) = x (mod 24). Usage: p309_census.py 10000000."""
import math
import sys
from collections import Counter
import numpy as np

N = int(sys.argv[1])
L = 6 * N
ph = np.arange(L + 1, dtype=np.int32)
sv = np.ones(L + 1, dtype=np.bool_)
sv[:2] = False
for i in range(2, math.isqrt(L) + 1):
    if sv[i]:
        sv[i * i::i] = False
for p in np.nonzero(sv)[0]:
    p = int(p)
    v = ph[p::p]
    v -= v // p
del sv
sig = np.zeros(N + 1, dtype=np.int64)
S = math.isqrt(N)
for d in range(1, S + 1):
    sig[d::d] += d
for k in range(1, S + 1):
    hi = N // k
    if hi <= S:
        break
    sig[k * (S + 1)::k] += np.arange(S + 1, hi + 1, dtype=np.int64)

x = np.arange(2, N + 1, dtype=np.int64)
sx = sig[2:]
assert sx.max() <= L
f = sig[ph[2:N + 1].astype(np.int64)] - ph[sx].astype(np.int64)
F = f - x
comp = ph[2:N + 1].astype(np.int64) != x - 1
print('solutions:', x[F == 0].tolist())
for lo, hi in ((10**3, 10**4), (10**4, 10**5), (10**5, 10**6), (10**6, N + 1)):
    m = comp & (np.abs(F) <= 100) & (x >= lo) & (x < hi)
    print('composite x in [%.0e, %.0e): %d near misses |F| <= 100' % (lo, hi, m.sum()))
c = Counter(F[comp & (x % 2 == 0) & (np.abs(F) <= 20)].tolist())
print('even composite x, count of each F in -20..20:', [(v, c.get(v, 0)) for v in range(-20, 21, 2)])
print('f even for', int((f[comp] % 2 == 0).sum()), 'of', int(comp.sum()), 'composite x')
print('x mod 24 : #(f = x mod 24) / #x')
for r in range(24):
    sel = comp & (x % 24 == r)
    print('  %2d: %d / %d' % (r, int((f[sel] % 24 == r).sum()), int(sel.sum())))
