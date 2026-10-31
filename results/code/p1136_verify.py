# Puzzle 1136: F(p) = d(p-1), G(p) = sigma(p-1). Casini: G prime => F prime. Converse: F prime => p = F_(2^n)
# (or p = 3), d(p-1) = F_n and G(p) = M_(F_n); so the converse holds unless F_n, F_(2^n) prime and M_(F_n) composite.
from sympy import divisor_count, divisor_sigma, isprime, primerange
both = []
for p in primerange(3, 2 * 10**6):
    F, G = divisor_count(p - 1), divisor_sigma(p - 1)
    assert isprime(F) == isprime(G), p           # Q1 holds for every prime below 2e6
    if isprime(F): both.append((p, F, G))
print("primes p < 2e6 with F(p), G(p) prime:", both)
assert both == [(3, 2, 3), (5, 3, 7), (17, 5, 31), (65537, 17, 131071)]
# the four cases are p = F_0, F_1, F_2, F_4 with G = M_2, M_3, M_5, M_17
Fermat = lambda k: 2**(2**k) + 1
for (p, F, G), (j, r) in zip(both, [(0, 2), (1, 3), (2, 5), (4, 17)]):
    assert p == Fermat(j) and F == r and G == 2**r - 1
# 257 = F_3: d(256) = 9 is not prime, sigma(256) = 511 = 7 * 73
assert divisor_count(256) == 9 and divisor_sigma(256) == 511 == 7 * 73
print("Q1 verified for p < 2e6; the converse can only fail at p = F_(2^n) with F_n prime, n >= 33 (F_5..F_32 composite)")
