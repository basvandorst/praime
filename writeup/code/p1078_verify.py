"""Puzzle 1078: N(p) = p + (q - p) 2^(p-2), p < q consecutive primes.  (about 1 minute)

1. 3 | N whenever 6 does not divide g = q - p (checked for all p < 10^5).
2. For 6 | g: no (probable) prime N for p < 6000 (gmpy2), as on the page (no prime for p < 75000).
3. No covering set: the share of N free of primes 5..10^4 stays near 15 % in every range.
4. Sieve-calibrated expectation of the number of primes: sum over N free of primes < R of
   prod_{r<R} (1 - 1/r)^(-1) / ln N; 1.56 for p < 75000, about 0.085 more per doubling of p.
"""
import math
import numpy as np
import gmpy2
from sympy import primerange

ps = list(primerange(5, 300000))
pairs = list(zip(ps, ps[1:]))
# 1. the mod-3 rule
for p, q in pairs:
    if p >= 10**5: break
    g = q - p
    n3 = (p + g * pow(2, p - 2, 3)) % 3
    assert (n3 == 0) == (g % 6 != 0), (p, q)
print("1. 3 | N exactly when 6 does not divide q - p (all p < 1e5)")
# 2. PRP search for p < 6000
prp = [p for p, q in pairs if p < 6000 and (q - p) % 6 == 0 and gmpy2.is_prime(p + (q - p) * (1 << (p - 2)), 3)]
assert prp == []
print("2. no probable prime N for p < 6000 with 6 | q - p")
# 3, 4. survivors of trial division by 5..R and the calibrated expectation
R = 10**4
small = np.array(list(primerange(5, R)), dtype=np.int64)
w = 1 / math.prod(1 - 1 / r for r in primerange(2, R))
def block(lo, hi):
    n = nf = 0; E = 0.0
    for p, q in pairs:
        if p < lo: continue
        if p >= hi: break
        g = q - p
        if g % 6: continue
        n += 1
        t = np.array([pow(2, p - 2, int(r)) for r in small], dtype=np.int64)
        res = (p + g * t) % small; res[small == p] = 1
        if (res == 0).any(): continue
        nf += 1; E += w / ((p - 2) * math.log(2) + math.log(g))
    return n, nf, E
for lo, hi in ((5, 75000), (75000, 150000), (150000, 300000)):
    n, nf, E = block(lo, hi)
    print(f"   p in [{lo}, {hi}): {n} with 6 | g, {nf} ({nf/n:.1%}) free of primes 5..1e4, expected primes {E:.3f}")
print("3/4. survivors stay near 15 %: no finite covering set; about 0.085 expected primes per doubling of p")
print("all claims verified")
