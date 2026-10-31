# Problem 76: check a set of distinct odd semiprimes whose reciprocals sum to 1 (writeup/prob-76.md).
# Run: python3 prob76_verify.py [DATAFILE]     (default ../data/prob76_solution.txt; a few seconds; needs sympy)
import os, sys
from fractions import Fraction
from sympy import factorint, primerange

path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), '..', 'data', 'prob76_solution.txt')
S = [int(x) for x in open(path).read().split()]

# 1. Distinct odd semiprimes (two prime factors, counted with multiplicity), and the exact sum.
assert len(S) == len(set(S))
fac = {}
for n in S:
    f = factorint(n)
    assert n % 2 == 1 and sum(f.values()) == 2, n
    assert len(f) == 2, n                      # squarefree: forced by the lemma of step 2
    fac[n] = sorted(f)
total = sum(Fraction(1, n) for n in S)
assert total == 1
print('%d distinct odd squarefree semiprimes, sum of reciprocals = %s, largest %d' % (len(S), total, max(S)))

# 2. The local criterion: for every prime p of the set, sum of q^-1 (mod p) over the partners q of p is 0.
nb = {}
for n in S:
    p, q = fac[n]
    nb.setdefault(p, []).append(q)
    nb.setdefault(q, []).append(p)
for p, qs in nb.items():
    assert sum(pow(q, -1, p) for q in qs) % p == 0
    assert len(qs) >= 2                        # degree 1 is impossible (q^-1 is never 0 mod p)
print('every one of the %d primes is balanced; minimum degree %d' % (len(nb), min(len(v) for v in nb.values())))

# 3. Why no p^2 can occur: with p^2 in the set the p-part of the sum is (p*A + 1)/p^2 with p not dividing p*A + 1.
#    Illustration: adding 9 = 3^2 to the solution leaves 3^2 in the denominator.
t = total + Fraction(1, 9)
assert t.denominator % 9 == 0

# 4. Lower bound: the 551 smallest odd squarefree semiprimes are the first whose reciprocals reach 1.
small, s, k = [], Fraction(0), 0
n = 3
while s < 1:
    f = factorint(n)
    if sum(f.values()) == 2 and len(f) == 2:
        s += Fraction(1, n)
        k += 1
    n += 2
assert k == 551
print('lower bound: any solution has at least %d terms (the smallest %d already reach 1 only at %d)' % (k, k, n - 2))
print('all checks passed')
