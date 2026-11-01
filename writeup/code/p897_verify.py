# Puzzle 897 Q2: check a set of distinct semiprimes whose reciprocals sum to an integer (writeup/897.md).
# Run: python3 p897_verify.py [DATAFILE] [TARGET]   (default ../data/p897_sum2.txt and 2; a few seconds)
import os, sys
from fractions import Fraction
import numpy as np
from sympy import factorint, primerange

here = os.path.dirname(__file__)
path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, '..', 'data', 'p897_sum2.txt')
target = int(sys.argv[2]) if len(sys.argv) > 2 else 2
S = [int(x) for x in open(path).read().split()]

# 1. Distinct semiprimes (squarefree: a square p^2 would leave p^2 in the denominator), exact sum.
assert len(S) == len(set(S))
for n in S:
    f = factorint(n)
    assert sum(f.values()) == 2 and len(f) == 2, n
total = sum(Fraction(1, n) for n in S)
assert total == target
print('%d distinct squarefree semiprimes; reciprocals sum to exactly %s; smallest %d, largest %d'
      % (len(S), total, min(S), max(S)))

# 2. Balance at every prime (equivalent to integrality).
nb = {}
for n in S:
    p, q = sorted(factorint(n))
    nb.setdefault(p, []).append(q); nb.setdefault(q, []).append(p)
assert all(sum(pow(q, -1, p) for q in qs) % p == 0 for p, qs in nb.items())
print('all %d primes are balanced' % len(nb))

# 3. Lower bound: the smallest semiprimes needed to reach the target (1234 for 2, 73736 for 3).
N = 400000
P = list(primerange(2, N // 2 + 1))
sp = sorted(p * q for i, p in enumerate(P) for q in P[i + 1:] if p * q <= N)
cs = np.cumsum(1.0 / np.array(sp, dtype=float))
k = int(np.searchsorted(cs, target)) + 1
assert sum(Fraction(1, x) for x in sp[:k]) >= target > sum(Fraction(1, x) for x in sp[:k - 1])
print('lower bound: at least %d terms' % k)
print('all checks passed')
