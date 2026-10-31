"""Puzzle 1181 Q2: F(p,q) = floor(q^2/p), q > p primes. Checks:
 (1) every prime F(p,q) (p < 600, p < q < 20000) has r = q^2 mod p even, and either r is not a square
     or q = p + s with r = s^2 (s even), F = p + 2s;
 (2) T(p) is empty for p = 2, 3, 5, 11 (q < 10^6; proved in writeup/1181.md);
 (3) every prime 7 <= p < 10^6 other than 11 has an even non-square quadratic residue r < p,
     and the case argument picks one: r = 2, or (p = 1 mod 4) counting, or r = 2(p - 4k^2);
 (4) T(p) is non-empty for every prime p < 3000 except 2, 3, 5, 11."""
import numpy as np, gmpy2
from math import isqrt
from sympy import primerange, isprime
def issq(n): return isqrt(n) ** 2 == n
for p in primerange(2, 600):
    for q in primerange(p + 1, 20000):
        F = q * q // p
        if isprime(F):
            r = q * q % p
            assert r % 2 == 0, (p, q)
            if issq(r):
                s = isqrt(r); assert q == p + s and F == p + 2 * s, (p, q)
for p in (2, 3, 5, 11):
    assert not any(gmpy2.is_prime(q * q // p) for q in primerange(p + 1, 10**6)), p
def legendre(a, p): return pow(a, (p - 1) // 2, p)
for p in primerange(7, 10**6):
    if p == 11: continue
    if p % 8 in (1, 7): r = 2                                      # 2 is a residue
    elif p % 4 == 1:                                               # p = 5 mod 8: half the residues are even
        r = next(x for x in range(2, p, 2) if legendre(x, p) == 1 and not issq(x))
    else:                                                          # p = 3 mod 8
        ks = [k for k in range(1, isqrt(p) + 1) if p / 2 < 4 * k * k < p]
        if ks: r = 2 * (p - 4 * ks[0] ** 2)                        # 2x, x odd non-residue < p/2
        else: r = next(x for x in range(2, p, 2) if legendre(x, p) == 1 and not issq(x))
    assert r < p and r % 2 == 0 and not issq(r) and legendre(r, p) == 1, (p, r)
small38 = [p for p in primerange(7, 10**6) if p % 8 == 3 and not any(p / 2 < 4*k*k < p for k in range(1, isqrt(p) + 1))]
print("p = 3 mod 8 with no k in (sqrt(p/8), sqrt(p/4)), handled directly:", small38)
for p in primerange(2, 3000):
    hit = next((q for q in primerange(p + 1, 10**6) if gmpy2.is_prime(q * q // p)), None)
    assert (hit is None) == (p in (2, 3, 5, 11)), p
print("all checks passed")
