# Puzzle 1178 (writeup/1178.md): F(n) = floor((p_n + 1)^2 / (p_{n+1} - 1)). With g = p_{n+1} - p_n,
#   (p + 1)^2 = (p - g + 3)(p + g - 1) + (g - 2)^2,   so   F(n) = p - g + 3 + floor((g - 2)^2 / (p + g - 1)),
# and p - g + 3 is even for p > 2. So F(n) is odd only if (g - 2)^2 >= p + g - 1, i.e. roughly g >= sqrt(p).
# Checks: the identity; F(n) prime only for n = 30 among all n with p_n < 2.2 * 10^6; the gap bound: every gap below
# 4 * 10^18 is at most 1476 (fact 3023 in FACTS.md, Oliveira e Silva-Herzog-Pardi 2014), and (g - 2)^2 < p + g - 1
# for g <= 1476 and p >= 2.2 * 10^6; the merit a further prime would need beyond 4 * 10^18.
# Run: python3 p1178_verify.py   (about 20 seconds; needs sympy, numpy)
from math import log, sqrt, isqrt
import numpy as np
from sympy import isprime, symbols, expand

p, g = symbols('p g')
assert expand((p + 1) ** 2 - ((p - g + 3) * (p + g - 1) + (g - 2) ** 2)) == 0
print('identity (p+1)^2 = (p-g+3)(p+g-1) + (g-2)^2 holds')

LIM = 2_200_000
s = np.ones(LIM + 1, dtype=bool); s[:2] = False
for i in range(2, isqrt(LIM) + 1):
    if s[i]: s[i * i::i] = False
P = np.nonzero(s)[0]
found = []
for n in range(len(P) - 1):
    pn, pn1 = int(P[n]), int(P[n + 1])
    F = (pn + 1) ** 2 // (pn1 - 1)
    gg = pn1 - pn
    assert F == pn - gg + 3 + (gg - 2) ** 2 // (pn + gg - 1)
    if isprime(F): found.append((n + 1, pn, pn1, F))
print('primes F(n) with p_n < %d: %s' % (LIM, found))
assert found == [(30, 113, 127, 103)]
GMAX = 1476                                   # largest gap below 4 * 10^18 (fact 3023)
assert all((gg - 2) ** 2 < 2_200_000 + gg - 1 for gg in range(2, GMAX + 1))
print('for every gap g <= %d and p >= 2.2e6: (g-2)^2 < p + g - 1, so F(n) is even; hence no other prime F(n) with p_n < 4e18' % GMAX)
x = 4e18
print('beyond 4e18 an odd F(n) needs g > sqrt(p) - 1, a gap of merit g/ln p > %.1e (the record merit is 42.5)' % ((sqrt(x) - 1) / log(x)))
print('all checks passed')
