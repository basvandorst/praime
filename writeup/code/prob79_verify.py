# Problem 79: check a set of numbers with exactly three prime factors whose reciprocals sum to 1 (writeup/prob-79.md).
# Run: python3 prob79_verify.py [DATAFILE]   (default ../data/prob79_solution.txt; about a second)
import os, sys
from fractions import Fraction
from sympy import factorint

path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), '..', 'data', 'prob79_solution.txt')
S = [int(x) for x in open(path).read().split()]

# 1. Distinct, exactly three prime factors counted with multiplicity, exact sum 1.
assert len(S) == len(set(S))
for n in S:
    assert sum(factorint(n).values()) == 3, n
assert sum(Fraction(1, n) for n in S) == 1
nsq = sum(1 for n in S if len(factorint(n)) == 2)
print('%d distinct numbers with three prime factors (%d not squarefree); reciprocals sum to exactly 1; largest %d'
      % (len(S), nsq, max(S)))

# 2. The criterion: R_p = p^2 * sum(1/t, p | t) is 0 mod p^2 for every prime p.
R = {}
for t in S:
    for p, e in factorint(t).items():
        c = pow(t // (p * p), -1, p * p) if e == 2 else (p * pow(t // p % p, -1, p)) % (p * p)
        assert e in (1, 2)                        # no cube p^3 in the set
        R[p] = (R.get(p, 0) + c) % (p * p)
assert all(r == 0 for r in R.values())
print('R_p = 0 (mod p^2) for all %d primes' % len(R))

# 3. Lower bound: the 189 smallest numbers with three prime factors (cubes excluded) are the first to reach 1.
cands = [n for n in range(8, 2000) if sum(factorint(n).values()) == 3 and len(factorint(n)) > 1]
s, k = Fraction(0), 0
for n in cands:
    s += Fraction(1, n); k += 1
    if s >= 1:
        break
assert k == 189 and n == 782
print('lower bound: at least %d terms (the %d smallest reach 1 at %d)' % (k, k, n))
print('all checks passed')
