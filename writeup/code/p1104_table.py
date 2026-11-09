# Puzzle 1104: least prime p for which S_n = sum of the n-th powers of the primitive roots of p is prime, for n = N0..N1 and p <= PMAX.
# Skips the p where S_n is provably composite: p | S_n when (p-1)/gcd(n, p-1) is not squarefree, and S_n is even when the number
# of odd primitive roots is even (both proved in the write-up). Primality: gmpy2.is_prime (PRP); the write-up's validation proves them.
# Usage: .venv/bin/python3 writeup/code/p1104_table.py N0 N1 PMAX
import sys, gmpy2
from sympy import primerange, factorint
from math import gcd
N0, N1, PMAX = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
ns = list(range(N0, N1 + 1)); best = {}
def squarefree(m): return all(e == 1 for e in factorint(m).values())
for p in primerange(3, PMAX + 1):
    todo = [n for n in ns if n not in best and squarefree((p - 1)//gcd(n, p - 1))]
    if not todo: continue
    m = p - 1; f = list(factorint(m)); g = 2
    while any(pow(g, m//q, p) == 1 for q in f): g += 1
    roots = []; x = 1
    for i in range(m):
        if all(i % q for q in f): roots.append(x)
        x = x*g % p
    if p > 3 and sum(r & 1 for r in roots) % 2 == 0: continue      # every S_n is even (and > 2)
    for n in todo:
        S = sum(pow(r, n) for r in roots)
        if gmpy2.is_prime(S): best[n] = (p, S); print(n, p, S, flush=True)
    if len(best) == len(ns): break
for n in ns:
    if n not in best: print(n, "none below", PMAX)
