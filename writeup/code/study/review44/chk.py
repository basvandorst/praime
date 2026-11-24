import numpy as np, math, gmpy2, random
SRC = open('mywin.py').read(); exec(SRC[SRC.index('def small_primes'):SRC.index('def block_stats')])
base = small_primes(10**6 + 10)
for lo in (6*10**9, 8399800000, 10**12 - 10**5):
    hi = lo + 2*10**5; P = sieve(lo, hi, base)
    G = [n for n in range(lo, hi) if gmpy2.is_prime(n)]
    print(lo, len(P), len(G), list(P) == G)
