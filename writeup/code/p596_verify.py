# Puzzle 596: why 12n(n+2)+1 is so often prime or semiprime (writeup/596.md).
# Run: python3 p596_verify.py   (about 10 seconds; needs sympy and numpy)
from math import log, sqrt
from collections import Counter
import numpy as np
from sympy import factorint, primerange, legendre_symbol, jacobi_symbol

f = lambda n: 12 * n * (n + 2) + 1

# 1. f(n) = 12m^2 - 11 with m = n + 1.  A prime q not dividing 66 divides some value iff (6m)^2 = 33 (mod q) is
#    solvable, i.e. (33/q) = 1, which by reciprocity is the Jacobi symbol (q/33) = 1.  2 and 3 never divide;
#    11 divides exactly when 11 | m, and then only once.
for n in range(0, 200):
    assert f(n) == 12 * (n + 1)**2 - 11
can = [q for q in primerange(2, 200) if any(f(n) % q == 0 for n in range(q))]
for q in primerange(5, 200):
    if q != 11:
        assert (q in can) == (legendre_symbol(33 % q, q) == 1) == (jacobi_symbol(q, 33) == 1)
assert 2 not in can and 3 not in can
assert [q for q in primerange(2, 60) if q not in can] == [2, 3, 5, 7, 13, 19, 23, 43, 47, 53, 59]
assert all((f(n) % 11 == 0) == ((n + 1) % 11 == 0) and f(n) % 121 != 0 for n in range(20000))

# 2. The smallest possible value with 3 prime factors is 11 * 17^2 = 3179 (11 only once; 17 is the least other
#    allowed prime), so f(n) is prime or semiprime for every n <= 15 (f(15) = 3061 < 3179 <= f(16) = 3457).
assert can[:2] == [11, 17] and f(15) < 11 * 17**2 <= f(16)
assert all(sum(factorint(f(n)).values()) <= 2 for n in range(16))

# 3. The page's data: exceptions (3 or more prime factors) for n <= 100, and per 1000 up to 10^4.
exc, om = [], Counter()
for n in range(0, 10000):
    k = sum(factorint(f(n)).values())
    om[(n // 1000, min(k, 3))] += 1
    if k >= 3:
        exc.append(n)
assert [n for n in exc if n <= 100] == [21, 38, 71, 76, 96]
assert f(21) == 5797 == 11 * 17 * 31
per = [sum(1 for n in exc if b * 1000 <= n < (b + 1) * 1000) for b in range(10)]
assert per == [180, 255, 279, 273, 309, 321, 304, 313, 333, 337]

# 4. The model: local factor omega(q) = number of roots of f mod q (0, 1 or 2); Hardy-Littlewood constant
#    C = prod (1 - omega(q)/q)/(1 - 1/q); P(prime) = C/ln f(n); P(semiprime) = sum over allowed p <= sqrt f of
#    (omega(p)/p) * (C / local factor at p) / ln(f/p).
def omega(q):
    return 0 if q in (2, 3) else 1 if q == 11 else (2 if jacobi_symbol(q, 33) == 1 else 0)
PR = list(primerange(2, 10**6))
C = 1.0
for q in PR:
    C *= (1 - omega(q) / q) / (1 - 1 / q)
assert abs(C - 4.81) < 0.01
A = np.array([q for q in PR if omega(q) > 0], dtype=float)
w = np.array([omega(int(q)) for q in A], dtype=float) / A
loc = (1 - w) / (1 - 1 / A)
def model(lo, hi, step=10):
    P = sum(C / log(f(n)) for n in range(max(lo, 1), hi))
    S = 0.0
    for n in range(max(lo, 1), hi, step):
        N = f(n)
        k = np.searchsorted(A, sqrt(N), side='right')
        S += float(np.sum(w[:k] * (C / loc[:k]) / np.log(N / A[:k]))) * step
    return P, S
print(' block of n         primes (model/data)   semiprimes (model/data)   3+ factors (model/data)')
for b in list(range(10)) + [100]:
    lo, hi = b * 1000, (b + 1) * 1000
    if b < 10:
        dp, ds = om[(b, 1)], om[(b, 2)]
    else:
        cnt = Counter(min(sum(factorint(f(n)).values()), 3) for n in range(lo, hi))
        dp, ds = cnt[1], cnt[2]
    P, S = model(lo, hi)
    E = (hi - max(lo, 1)) - P - S
    de = (hi - lo) - dp - ds
    print(' [%6d, %6d)   %6.1f / %4d          %6.1f / %4d              %6.1f / %4d' % (lo, hi, P, dp, S, ds, E, de))
    assert abs(P - dp) < 40 and abs(S - ds) < 40 and abs(E - de) < 40     # within about 2.5 standard deviations
print('all checks passed')
