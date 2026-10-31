"""Puzzle 24: checks for writeup/24.md.
 (1) the four known primes that stay prime read in bases 2..9 (BPSW: probable primes above 2^64) and Resta's base-2..11 prime;
 (2) the model against an exact count: primes n in [1e7, 3e7] with M_7, M_8, M_9 all prime;
 (3) the Monte Carlo model below 5e10, 2.8e11 and 3.1e12 (bases 2..9);
 (4) the analytic model and its calibration factor (Monte Carlo / analytic), then the base-2..12 forecast.
Runtime about 25 minutes (the Monte Carlo parts). Run from writeup/code."""
import math, gmpy2
from sympy import isprime
import p24_model_mc as mc, p24_model_analytic as an
def Mb(n, b):
    s, p = 0, 1
    while n: n, d = divmod(n, b); s += d * p; p *= 10
    return s
for n in (50006393431, 727533146383, 2250332130313, 2651541199513):
    assert isprime(n) and all(isprime(Mb(n, b)) for b in range(2, 10)), n
n = 181395559296673
assert isprime(n) and all(isprime(Mb(n, b)) for b in range(2, 12)) and all(d < 10 for d in [n // 11**i % 11 for i in range(15)])
P = mc.primes_window(10**7, 2 * 10**7)
real = sum(1 for n in P.tolist() if all(gmpy2.is_prime(Mb(n, b)) for b in (7, 8, 9)))
print(f"exact check, bases 7,8,9 on [1e7,3e7]: real {real}, model {mc.count(1e7, 3e7, [7, 8, 9], npts=3):.1f}")
B9 = list(range(2, 10))
mcN = {X: mc.count(1e8, X, B9) for X in (5.0006e10, 2.8e11, 3.1e12)}
for X, N in mcN.items():
    print(f"bases 2..9: expected below {X:.3g}: {N:.2f}  P(at least one) = {1 - math.exp(-N):.2f}")
cal = mcN[3.1e12] / an.count(1e6, 3.1e12, list(range(2, 11)))[0]
print(f"calibration factor Monte Carlo / analytic: {cal:.2f}")
B12 = list(range(2, 13))
for X in (6e15, 1e17, 1e18, 1e19):
    N = cal * an.count(1e6, X, B12)[0]
    print(f"bases 2..12: expected below {X:.3g}: {N:.3f}  P(at least one) = {1 - math.exp(-N):.2f}")
