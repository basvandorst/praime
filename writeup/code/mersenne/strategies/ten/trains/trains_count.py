"""Exact number of trains (runs of consecutive primes) with sum <= x, for x up to 10^9: for every start i the longest run
is found by a binary search in the prefix sums (numpy searchsorted); primes from a segmented sieve.  The mean of f(n)
over n <= x is that number divided by x.   usage: trains_count.py out.txt"""
import sys, math
import numpy as np
out = open(sys.argv[1], 'w')
def pr(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')
def primes_upto(X):
    r = int(X ** .5) + 1; s = np.ones(r + 1, bool); s[:2] = False
    for i in range(2, int(r ** .5) + 1):
        if s[i]: s[i * i::i] = False
    base = np.nonzero(s)[0]; parts = [base]; seg = 10 ** 7
    for lo in range(r + 1, X + 1, seg):
        hi = min(X + 1, lo + seg); m = np.ones(hi - lo, bool)
        for p in base:
            st = max(p * p, -(-lo // p) * p)
            if st >= hi: continue
            m[st - lo::p] = False
        parts.append((np.nonzero(m)[0] + lo).astype(np.int64))
    return np.concatenate(parts).astype(np.int64)
pr('        x   trains with sum <= x   mean f      (ln 2 = %.6f)' % math.log(2))
P = primes_upto(10 ** 9)
cs = np.concatenate([[0], np.cumsum(P)])
for e in range(3, 10):
    x = 10 ** e; n = int(np.searchsorted(P, x, 'right')); c = cs[:n + 1]
    j = np.searchsorted(c, c[:n] + x, 'right') - 1          # largest j with c[j] - c[i] <= x
    cnt = int(np.sum(j - np.arange(n)))
    pr('   10^%d %18d %10.6f' % (e, cnt, cnt / x))
