"""Puzzle 1166, independent check.

1. Direct product: for each (n, X) on the page and in our table, prod_{x in X} P_n(x) is computed
   literally (product of p (p + x) over pairs of odd primes p < p + x <= p_n) and compared with
   (p_n# / 2)^4 as plain integers.
2. Independent solver: for n = 3 .. NMAX the 0/1 system 'every odd prime q <= p_n has exactly
   four partners q +- x among the odd primes <= p_n' is given to HiGHS (scipy.optimize.milp).
   Feasible n must be exactly {6, 8, 9, 11, 13, 14, 15}; infeasible otherwise.
Usage: p1166_validate.py [NMAX]
"""
import sys
import numpy as np
from sympy import prime, primerange, primorial
from scipy.optimize import milp, LinearConstraint, Bounds

page = {6: [2, 4, 6, 8, 10], 8: [2, 4, 6, 14, 16], 9: [4, 8, 10, 12, 14, 18],
        11: [6, 8, 14, 16, 22, 24, 28], 13: [6, 14, 16, 22, 26, 28, 30, 36],
        14: [6, 12, 16, 28, 32, 34, 38], 15: [6, 12, 16, 32, 34, 38, 44]}

def P(n, x):
    pn = prime(n); odd = list(primerange(3, pn + 1)); S = set(odd); r = 1
    for p in odd:
        if p + x in S: r *= p * (p + x)
    return r

for n, X in page.items():
    prod = 1
    for x in X: prod *= P(n, x)
    assert prod == (primorial(n) // 2) ** 4, n
    print(f"n = {n}: prod P_n(x) over {X} = (p_n#/2)^4  ok")

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 60
feasible = []
for n in range(3, NMAX + 1):
    pn = prime(n); odd = list(primerange(3, pn + 1)); S = set(odd)
    xs = list(range(2, pn - 2, 2))
    A = np.array([[(q - x in S) + (q + x in S) for x in xs] for q in odd], dtype=float)
    res = milp(c=np.zeros(len(xs)), constraints=LinearConstraint(A, 4, 4),
               integrality=np.ones(len(xs)), bounds=Bounds(0, 1))
    if res.status == 0:
        X = [x for x, v in zip(xs, res.x) if v > 0.5]
        feasible.append(n)
        assert n in page and X == page[n], (n, X)
    else:
        assert res.status == 2, (n, res.status, res.message)   # 2 = infeasible
print(f"HiGHS, n = 3..{NMAX}: feasible exactly for n in {feasible}")
assert feasible == sorted(page)
