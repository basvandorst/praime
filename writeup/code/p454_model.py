"""Puzzle 454: Hardy-Littlewood H(x) = (C/2)(Li(sqrt x) - Li(sqrt 2)) against the exact counts
(e <= 16 from p454_count.py 100000000; e = 17..28 from OEIS A083844, not recomputed), the page's three approximations, and the crossing x*."""
from mpmath import mp, mpf, li, sqrt, log, exp
import numpy as np

mp.dps = 30
C = mpf('1.3728134628182460091')
exact = {1: 2, 2: 4, 3: 10, 4: 19, 5: 51, 6: 112, 7: 316, 8: 841, 9: 2378, 10: 6656, 11: 18822,
         12: 54110, 13: 156081, 14: 456362, 15: 1339875, 16: 3954181}
oeis = {17: 11726896, 18: 34900213, 19: 104248948, 20: 312357934, 21: 938457801, 22: 2826683630,
        23: 8533327397, 24: 25814570672, 25: 78239402726, 26: 237542444180, 27: 722354138859,
        28: 2199894223892}
H = lambda x: C / 2 * (li(sqrt(x)) - li(sqrt(2)))
G = lambda x: C * sqrt(x) / log(x)
fa = lambda x: sqrt(2 * x) / log(x)
fb = lambda x: sqrt(2 * x) / (log(x) - 1)
fc = lambda x: sqrt(2 / x) * li(x)

ps = np.arange(3, 10**7, 2)
sv = np.ones(10**7, dtype=bool); sv[:2] = False
for i in range(2, 3163):
    if sv[i]: sv[i * i::i] = False
ps = np.nonzero(sv)[0][1:].astype(np.float64)
chi = np.where(ps % 4 == 1, 1.0, -1.0)
acc = np.where(chi > 0, 1 - 1 / (ps - 1) ** 2, 1 + 1 / (ps ** 2 - 1))
print('C = (4/pi) * absolutely convergent product over odd p < 1e7: %.9f (quoted %.9f)' % (4 / np.pi * np.prod(acc), C))

print(' e      exact   (H-exact)/sqrt(exact)   G/ex-1   a/ex-1   b/ex-1   c/ex-1')
for e, n in {**exact, **oeis}.items():
    x = mpf(10) ** e
    print('%2d %10d %14.2f %16.4f %8.4f %8.4f %8.4f' % (e, n, (H(x) - n) / sqrt(n), G(x) / n - 1,
                                                       fa(x) / n - 1, fb(x) / n - 1, fc(x) / n - 1))

def bisect(f, lo, hi):
    for _ in range(200):
        mid = (lo + hi) / 2
        if f(mid) * f(lo) > 0:
            lo = mid
        else:
            hi = mid
    return lo

for name, f in (('a', fa), ('b', fb), ('c', fc)):
    d = lambda L: abs(f(exp(L)) / H(exp(L)) - 1) - abs(G(exp(L)) / H(exp(L)) - 1)
    print(name, 'is closer to H than C sqrt(x)/ln x from 1e20 up to x = %.3e' % exp(bisect(d, mpf(46), mpf(400))))
Ls = bisect(lambda L: fa(exp(L)) - H(exp(L)), mpf(60), mpf(80))
print('x* (sqrt(2x)/ln x = H) = %.6e' % exp(Ls))
for k in (1, 10):
    lo = bisect(lambda L: H(exp(L)) - fa(exp(L)) - k * sqrt(H(exp(L))), mpf(46), Ls)
    hi = bisect(lambda L: fa(exp(L)) - H(exp(L)) - k * sqrt(H(exp(L))), Ls, mpf(80))
    print('|a - H| < %d sqrt(H) only on [%.5e, %.5e]' % (k, exp(lo), exp(hi)))
x = mpf(10) ** 20
print('at 1e20: H - a = %.3e = %.0f sqrt(H)' % (H(x) - fa(x), (H(x) - fa(x)) / sqrt(H(x))))
r = [fa(exp(mpf(46) + i / mpf(4))) / H(exp(mpf(46) + i / mpf(4))) for i in range(2001)]
print('a/H increasing for ln x in [46, 546]:', all(u < v for u, v in zip(r, r[1:])), '; limit sqrt(2)/C = %.6f' % (sqrt(2) / C))
