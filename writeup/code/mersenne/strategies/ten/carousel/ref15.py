"""Independent reference for prime15.html: for each base a, the number of primes p <= X with p not dividing a, and how
many of them have a as a primitive root (sympy.is_primitive_root).   usage: ref15.py X bases out.json"""
import sys, json
from sympy import primerange, is_primitive_root
X = int(sys.argv[1]); bases = [int(b) for b in sys.argv[2].split(',')]; out = {}
P = list(primerange(2, X + 1))
for a in bases:
    t = c = 0
    for p in P:
        if a % p == 0: continue
        t += 1; c += is_primitive_root(a % p, p)
    out[str(a)] = [c, t]
json.dump(out, open(sys.argv[3], 'w')); print(out)
