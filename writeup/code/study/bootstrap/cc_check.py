# Compare with the published numbers of Connes-Consani (arXiv 2106.01715): smallest eigenvalue at mu = 3
# ("< 6e-8", even part), at lambda^2 = 11 ("2.389e-48"), and positivity lost at p = 1.9999 and 2.0005 (mu = 3).
import sys, time
from mpmath import mp, mpf, log, eigsy, nstr
from sympy import primerange
from weil import build_A, comb, sectors
def low(X, N, pos):
    a = log(mpf(X))/2
    E, O = sectors(a, N, build_A(a, N) - comb(a, N, pos))
    return min(eigsy(E, eigvals_only=True)), min(eigsy(O, eigvals_only=True))
mp.dps = 40
for N in (12, 24, 48):
    print('mu=3 N=%d even %s odd %s' % (N, *[nstr(v, 5) for v in low(3, N, [2])]), flush=True)
for p in ('1.9999', '2.0005'):
    print('mu=3 one point at', p, ': lam_min', nstr(min(low(3, 24, [mpf(p)])), 4), flush=True)
mp.dps = 80
for N in (24, 36, 48):
    t = time.time()
    e, o = low(11, N, [p for p in primerange(2, 11)])
    print('mu=11 N=%d even %s odd %s (%.0fs)' % (N, nstr(e, 5), nstr(o, 5), time.time() - t), flush=True)
