"""Independent reference for prime21.html (sympy's legendre_symbol and primitive_root): for the odd primes p <= X, the
count of each least quadratic non-residue, the sum of the least non-residues and of the least primitive roots.
usage: ref21.py X out.json"""
import sys, json
from sympy import primerange, legendre_symbol, primitive_root
X = int(sys.argv[1]); cnt = {}; s = g = 0; n = 0
for p in primerange(3, X + 1):
    q = 2
    while legendre_symbol(q, p) == 1: q += 1
    cnt[q] = cnt.get(q, 0) + 1; s += q; g += primitive_root(p); n += 1
json.dump(dict(n=n, cnt={str(k): v for k, v in sorted(cnt.items())}, sum_n=s, sum_g=g), open(sys.argv[2], 'w')); print(n, s / n, g / n)
