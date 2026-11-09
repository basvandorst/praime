"""Puzzle 58 / 54, Python version of p58_search.c for small ranges: consecutive primes p in [LO, HI) with
4p^2 + 1 prime (exact, gmpy2), report runs of >= K.  Usage: p58_search.py LO HI K"""
import sys
import numpy as np
from gmpy2 import mpz, is_prime
LO, HI, K = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
# segmented sieve of [LO, HI)
R = int(HI ** 0.5) + 1
s = np.ones(R + 1, bool); s[:2] = False
for i in range(2, int(R ** 0.5) + 1):
    if s[i]: s[i*i::i] = False
base = np.flatnonzero(s)
seg = np.ones(HI - LO, bool)
for p in base:
    st = max(p * p, ((LO + p - 1) // p) * p); seg[st - LO::p] = False
primes = np.flatnonzero(seg) + LO
run = []
for p in primes.tolist():
    if is_prime(mpz(4 * p * p + 1)): run.append(p)
    else:
        if len(run) >= K: print(f"k={len(run)} run starting at {run[0]}: {run}")
        run = []
print(f"{len(primes)} primes in [{LO}, {HI}) checked")
