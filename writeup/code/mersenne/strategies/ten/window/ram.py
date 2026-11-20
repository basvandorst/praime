"""Ramanujan primes: R_n is the least number with pi(x) - pi(x/2) >= n for all x >= R_n (Ramanujan 1919, Sondow 2009).
Share of the primes up to x that are Ramanujan primes, against the limit 1/2 (Amersi, Beckwith, Miller, Ronan and Sondow
2013: the share tends to 1/2), and R_n against p_{2n}.  usage: ram.py X out.txt"""
import sys, math
import numpy as np
X = int(sys.argv[1]); out = open(sys.argv[2], 'w')
def pr(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')
sv = np.ones(X + 1, bool); sv[:2] = False
for i in range(2, int(X ** .5) + 1):
    if sv[i]: sv[i * i::i] = False
pi = np.cumsum(sv)
x = np.arange(X + 1); s = pi - pi[x // 2]                    # s(x) = pi(x) - pi(x/2)
# R_n = 1 + max{x <= X : s(x) <= n - 1}; valid while that max is far below X (s grows), use n <= s(X)/2 for safety
nmax = int(s[X] // 2)
last = np.zeros(nmax + 2, np.int64)                            # last[k] = largest x with s(x) <= k  (s(x) <= k  <=> ...)
mx = np.maximum.accumulate(s[::-1])[::-1]                      # mx[x] = max_{y >= x} s(y)
# largest x with s(x) = k for each k: scan
for k in range(0, nmax + 1):
    pass
idx = np.nonzero(s <= nmax)[0]
lastx = np.zeros(nmax + 1, np.int64)
order = np.argsort(s[idx], kind='stable'); vals = s[idx][order]; xs = idx[order]
for k in range(nmax + 1):
    lo = np.searchsorted(vals, k, 'right')
    lastx[k] = xs[:lo].max() if lo > 0 else 0
R = np.maximum.accumulate(lastx)[:-0 or None] + 1             # R_{k+1} = 1 + max{x : s(x) <= k}
R = R[:nmax]                                                   # R[n-1] = R_n, n = 1..nmax
isR = np.zeros(X + 1, bool); isR[R] = True
P = np.nonzero(sv)[0]
pr('X = %d: %d Ramanujan primes computed (R_1 = %d, R_2 = %d, R_3 = %d, R_%d = %d); all prime: %s' % (X, nmax, R[0], R[1], R[2], nmax, R[-1], bool(np.all(sv[R]))))
lim = R[-1]
for e in range(3, int(math.log10(lim)) + 1):
    y = 10 ** e
    if y > lim: break
    pr('  up to 10^%d: share of primes that are Ramanujan primes %.4f' % (e, isR[:y + 1].sum() / sv[:y + 1].sum()))
pr('  R_n / p_{2n} for n = 10, 100, 1000, 10000: ' + ', '.join('%.4f' % (R[n - 1] / P[2 * n - 1]) for n in (10, 100, 1000, 10000) if n <= nmax and 2 * n <= len(P)))
