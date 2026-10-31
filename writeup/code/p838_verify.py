"""Puzzle 838: prime Kaprekar numbers other than repunit primes.
k is Kaprekar iff k^2 = a*10^n + b with 0 < b < 10^n and a + b = k (Iannucci: k = 0 mod d, k = 1 mod N/d,
d a unitary divisor of N = 10^n - 1). For prime k this means: k = p | N with m = N/p dividing p - 1.
Checks: (1) brute force, Kaprekar condition vs that criterion for k < 2*10^5; (2) even n impossible (the proof's
two facts, n <= 60); (3) for odd n, p must be primitive (p | Phi_n(10)), so phi(n)*ln 11 >= (n/2) ln 10 is needed;
(4) search n <= 600 with trial division to 2e6: only the repunit R19 (m = 9) appears."""
import gmpy2
from sympy import isprime, primerange, totient, cyclotomic_poly
def kaprekar(k):
    s = k * k
    for n in range(1, len(str(s)) + 1):
        a, b = divmod(s, 10**n)
        if b > 0 and a + b == k: return True
    return False
for k in range(2, 2 * 10**5):
    if isprime(k):
        crit = any((10**n - 1) % k == 0 and k > (10**n - 1) // k and (k - 1) % ((10**n - 1) // k) == 0
                   and gmpy2.gcd(k, (10**n - 1) // k) == 1 for n in range(1, 2 * len(str(k)) + 1))
        assert kaprekar(k) == crit, k
        assert not kaprekar(k), k               # no prime Kaprekar number below 2e5
for k in range(1, 60):                          # even n = 2k: the only candidate p = 10^k + 1, m = 10^k - 1
    assert gmpy2.gcd(10**k - 1, 10**k) == 1 and 10**k - 1 > 1
odd_excluded = [n for n in range(3, 2000, 2) if totient(n) * gmpy2.log(11) < n / 2 * gmpy2.log(10)]
print("odd n <= 2000 excluded by phi(n)/n < 0.48:", len(odd_excluded), "e.g.", odd_excluded[:6])
SP = list(primerange(2, 2 * 10**6)); found = []
for n in range(1, 601):
    N = 10**n - 1; c = N; m = 1
    for q in SP:
        if q * q > c: break
        while c % q == 0: c //= q; m *= q
    if c > 1 and m < c and gmpy2.is_prime(c) and (c - 1) % m == 0: found.append((n, m))
assert found == [(19, 9)], found
print("n <= 600, smooth part up to 2e6: only n = 19 (the repunit R19); all checks passed")
