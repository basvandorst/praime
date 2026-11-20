"""Ramanujan primes: check of the identity pi_R(x) = min_{y >= x} (pi(y) - pi(y/2)) (counted over y up to 2*10^7), and
of the share pi_R(x) / pi(x) against the smooth value (li(x) - li(x/2)) / li(x) and against 1/2.  usage: ram2.py out.txt"""
import sys, math
import numpy as np
from scipy.special import expi
out = open(sys.argv[1], 'w')
def pr(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')
X = 20000000
sv = np.ones(X + 1, bool); sv[:2] = False
for i in range(2, int(X ** .5) + 1):
    if sv[i]: sv[i * i::i] = False
pi = np.cumsum(sv); x = np.arange(X + 1); s = pi - pi[x // 2]
fmin = np.minimum.accumulate(s[::-1])[::-1]                 # fmin[x] = min_{x <= y <= X} s(y)
# Ramanujan primes from the definition: R_n = 1 + max{y : s(y) <= n - 1}
nmax = int(fmin[X // 2]); last = np.full(nmax + 1, 0, np.int64)
order = np.argsort(s, kind='stable'); sv_s = s[order]
for k in range(nmax + 1):
    hi = np.searchsorted(sv_s, k, 'right'); last[k] = order[:hi].max() if hi else 0
R = np.maximum.accumulate(last)[:nmax] + 1
pr('%d Ramanujan primes from the definition, largest %d; all prime: %s' % (nmax, R[-1], bool(np.all(sv[R]))))
li = lambda t: expi(math.log(t)) - expi(math.log(2))
pr('       x    pi_R(x) counted   min_{y>=x} s(y)   share   (li(x)-li(x/2))/li(x)   1/2 - share   times ln x')
for e in (3, 4, 5, 6, 7):
    y = 10 ** e
    if y > R[-1]: break
    c = int(np.searchsorted(R, y, 'right')); m = int(fmin[y]); sh = c / pi[y]; sm = (li(y) - li(y / 2)) / li(y)
    pr('%8d %16d %17d %8.4f %20.4f %13.4f %11.3f' % (y, c, m, sh, sm, 0.5 - sh, (0.5 - sh) * math.log(y)))
