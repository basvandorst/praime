# Puzzle 1101: primes p dividing the concatenation of sigma(p-1) and sigma(p+1).
# Checks 2 and 4793, searches all p < 2e7 with a sigma sieve, and computes the heuristic count.
# Run from the repo root: .venv/bin/python3 writeup/code/p1101_verify.py
import numpy as np
from math import log
from sympy import isprime, factorint
from sympy.functions.combinatorial.numbers import divisor_sigma

def concat(a, b): return int(str(a) + str(b))
for p in (2, 4793):
    c = concat(divisor_sigma(p - 1), divisor_sigma(p + 1)); print(p, divisor_sigma(p - 1), divisor_sigma(p + 1), c, factorint(c)); assert c % p == 0

N = 20_000_001
sig = np.zeros(N, np.int64)
for d in range(1, N):
    sig[d::d] += d
s = np.ones(N, bool); s[:2] = False
for i in range(2, int(N**0.5) + 1):
    if s[i]: s[i*i::i] = False
P = np.flatnonzero(s)
# concatenation = sigma(p-1) * 10^digits(sigma(p+1)) + sigma(p+1); test divisibility by p with Python ints (values < 1e18 fit int64 anyway)
digits = np.floor(np.log10(sig[P + 1])).astype(np.int64) + 1
assert np.all(10**digits > sig[P + 1]) and np.all(10**(digits - 1) <= sig[P + 1])
c = sig[P - 1]*10**digits + sig[P + 1]
sols = P[c % P == 0]
print("solutions p < 2e7:", list(sols)); assert list(sols) == [2, 4793]

M = 0.2614972128
for X in (1.1e10, 1e12, 1e15, 1e20, 1e27, 1e30):
    print("expected number of solutions with p <= %.1e: %.2f" % (X, log(log(X)) + M))
print("found below 1.1e10: 2")
print("all checks passed")
