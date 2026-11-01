# Puzzle 1020: check a set of squarefree numbers with exactly three prime factors whose reciprocals sum to 1
# (writeup/1020.md).  Run: python3 p1020_verify.py [DATAFILE]   (default ../data/p1020_solution.txt; ~10 s)
import os, sys
from fractions import Fraction
import numpy as np
from sympy import factorint, primerange

path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), '..', 'data', 'p1020_solution.txt')
S = [int(x) for x in open(path).read().split()]

# 1. Distinct, squarefree, exactly three prime factors; exact sum 1.
assert len(S) == len(set(S))
fac = {}
for n in S:
    f = factorint(n)
    assert len(f) == 3 and all(e == 1 for e in f.values()), n
    fac[n] = sorted(f)
assert sum(Fraction(1, n) for n in S) == 1
print('%d distinct squarefree numbers with three prime factors; reciprocals sum to exactly 1; '
      'smallest %d, largest %d' % (len(S), min(S), max(S)))

# 2. The balance criterion holds at every prime (equivalent to the sum being an integer).
res = {}
for n in S:
    a, b, c = fac[n]
    for x, y, z in ((a, b, c), (b, a, c), (c, a, b)):
        res[x] = (res.get(x, 0) + pow(y * z % x, -1, x)) % x
assert all(r == 0 for r in res.values())
print('all %d primes are balanced' % len(res))

# 3. Lower bound: the 4402 smallest squarefree numbers with three prime factors are the first to reach 1.
N = 30000
omega = np.zeros(N + 1, np.int8)
sqf = np.ones(N + 1, bool)
for p in primerange(2, N + 1):
    omega[p::p] += 1
    if p * p <= N:
        sqf[p * p::p * p] = False
idx = np.nonzero((omega == 3) & sqf)[0]
s, k = Fraction(0), 0
for n in idx:
    s += Fraction(1, int(n)); k += 1
    if s >= 1:
        break
assert k == 4402 and int(n) == 23309
print('lower bound: at least %d terms (the %d smallest reach 1 at %d)' % (k, k, int(n)))
print('all checks passed')
