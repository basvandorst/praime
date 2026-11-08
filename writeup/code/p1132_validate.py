"""Puzzle 1132: check the chains for k = 5 and k = 6 and the upper bound that proves them longest.

Chain check: every entry is a k-digit prime (sympy.isprime is a proof below 2^64), all entries are
distinct, and neighbours differ in exactly one digit position.
Bound: build the graph of k-digit primes (edges = one digit changed) from a numpy sieve; an isolated
prime is in no chain of length > 1, and a prime of degree 1 can only be an end of a chain; so
K <= (number of primes) - (isolated) - (degree-1 primes - 2).
Run: python3 p1132_validate.py   (about 10 s)
"""
import os
import numpy as np
from sympy import isprime

here = os.path.dirname(os.path.abspath(__file__))

def sieve(N):
    s = np.ones(N, dtype=bool); s[:2] = False
    for i in range(2, int(N**0.5) + 1):
        if s[i]: s[i*i::i] = False
    return s

for k in (5, 6):
    chain = [int(x) for x in open(os.path.join(here, '..', 'data', f'p1132_k{k}.txt'))]
    assert all(10**(k-1) <= p < 10**k and isprime(p) for p in chain)
    assert len(set(chain)) == len(chain)
    for a, b in zip(chain, chain[1:]):
        sa, sb = str(a), str(b)
        assert sum(x != y for x, y in zip(sa, sb)) == 1, (a, b)
    s = sieve(10**k); primes = [p for p in range(10**(k-1), 10**k) if s[p]]
    deg = {}
    for p in primes:
        d, t = 0, str(p)
        for i in range(k):
            for c in '0123456789':
                if c != t[i] and not (i == 0 and c == '0') and s[int(t[:i] + c + t[i+1:])]: d += 1
        deg[p] = d
    iso = sum(1 for p in primes if deg[p] == 0); leaf = [p for p in primes if deg[p] == 1]
    bound = len(primes) - iso - (len(leaf) - 2)
    assert len(chain) == bound and deg[chain[0]] == 1 and deg[chain[-1]] == 1
    print(f"k = {k}: chain of {len(chain)} distinct primes, one digit changed per step; "
          f"{len(primes)} primes, {iso} isolated, {len(leaf)} of degree 1 -> bound {bound}: optimal")
    print(f"   ends {chain[0]} ... {chain[-1]}; isolated: {[p for p in primes if deg[p] == 0]}")
print("all checks passed")
