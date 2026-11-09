"""Puzzle 1049.  f(p) = sum of d^d over the digits d of p (zero-free primes).  chain(p) = number of
steps p -> f(p) -> ... while the values stay zero-free primes.  The page asks for the least zero-free
prime with chain >= n.  f depends only on the digit multiset, so we enumerate multisets of 1..9 by
size, compute the chain of f(multiset), and for the qualifying multisets the least prime arrangement.
Usage: p1049_search.py maxdigits n"""
import sys, itertools
from gmpy2 import mpz, is_prime
from sympy.utilities.iterables import multiset_permutations
POW = {d: d**d for d in range(1, 10)}
def f(x): return sum(POW[int(c)] for c in str(x))
def zerofree_prime(x): return '0' not in str(x) and is_prime(x)
def chain(p):
    """number of successive zero-free primes after p"""
    k = 0
    while True:
        q = f(p)
        if not zerofree_prime(q): return k
        k += 1; p = q
        if k > 20: return k
maxd = int(sys.argv[1]); n = int(sys.argv[2])
best = None
for size in range(1, maxd + 1):
    cands = []
    for ms in itertools.combinations_with_replacement(range(1, 10), size):
        S = sum(POW[d] for d in ms)
        if not zerofree_prime(S): continue
        if 1 + chain(S) >= n: cands.append(ms)
    found = []
    for ms in cands:
        # least prime arrangement: permutations in increasing order
        digs = sorted(ms)
        for perm in multiset_permutations(digs):
            v = int("".join(map(str, perm)))
            if is_prime(v): found.append((v, ms)); break
    if found:
        found.sort(); v, ms = found[0]
        print(f"n={n}: least prime with {size} digits: {v} (multiset {''.join(map(str, ms))}); {len(cands)} qualifying multisets of this size")
        p = v; seq = [p]
        for _ in range(n): p = f(p); seq.append(p)
        print(" chain:", " -> ".join(map(str, seq)))
        break
    print(f"size {size}: {len(cands)} qualifying multisets, no prime arrangement" if cands else f"size {size}: none")
