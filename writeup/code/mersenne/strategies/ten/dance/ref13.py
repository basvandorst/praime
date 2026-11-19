"""Independent reference for prime13.html: for q in 3, 4, 5, 8, 10, 12, the transition counts between the remainders
mod q of consecutive primes (both > q), up to X (plain Python, sympy primes).   usage: ref13.py X out.json"""
import sys, json
from math import gcd
from sympy import primerange
X = int(sys.argv[1]); P = list(primerange(2, X + 1)); out = {}
for q in (3, 4, 5, 8, 10, 12):
    spots = [r for r in range(1, q) if gcd(r, q) == 1]; ix = {r: i for i, r in enumerate(spots)}
    M = [[0] * len(spots) for _ in spots]; prev = None
    for p in P:
        if p <= q: continue
        j = ix[p % q]
        if prev is not None: M[prev][j] += 1
        prev = j
    tot = sum(map(sum, M)); out[str(q)] = dict(M=M, tot=tot, rep=round(sum(M[i][i] for i in range(len(M))) / tot, 6))
json.dump(out, open(sys.argv[2], 'w')); print({q: (v['tot'], v['rep']) for q, v in out.items()})
