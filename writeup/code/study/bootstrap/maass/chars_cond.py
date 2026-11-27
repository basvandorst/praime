# Exploration (not registered): real characters with the values at 2 and 3 prescribed.  For each parity and each
# (x_2, x_3) in {-1, 0, 1}^2: the smallest admissible conductor (one measure on [-1, 1] at the other primes) against the
# smallest conductor of a real primitive character of that parity with those values.  Usage: chars_cond.py <X>
import sys, numpy as np
from sympy import kronecker_symbol, isprime
from sympy.ntheory import factorint
from maass import Problem
X = float(sys.argv[1])
def fundamental(D):
    if D % 4 == 1: return all(e == 1 for e in factorint(abs(D)).values())
    if D % 4 == 0:
        m = D//4; return m % 4 in (2, 3) and all(e == 1 for e in factorint(abs(m)).values())
    return False
grid = np.unique(np.concatenate([np.linspace(-1, 1, 81), np.linspace(-1, -0.98, 21), np.linspace(0.98, 1, 21), np.linspace(-0.02, 0.02, 21)]))
for kappa in (1, 0):
    for x2 in (-1, 0, 1):
        for x3 in (-1, 0, 1):
            true = next(abs(D) for n in range(3, 2000) for D in ((-n,) if kappa else (n,)) if fundamental(D) and kronecker_symbol(D, 2) == x2 and kronecker_symbol(D, 3) == x3)
            cf = lambda ph: ph; cf.pm = lambda m, x: x**m
            P = Problem(0.0, 0, X, T=[('chi', [kappa], 1, cf)], th=grid, margin=14.0)
            mask = np.ones((len(P.primes), len(grid)))
            for i, v in ((0, x2), (1, x3)): mask[i] = (np.abs(grid - v) < 1e-12)
            lam, t, rd = P.solve(mask=mask, tol=1e-7, rounds=120)
            D0 = -true if kappa else true
            match = sum(abs(float(w @ P.th) - int(kronecker_symbol(D0, p))) < 0.02 for p, w in zip(P.primes, P.w))
            print('%s x_2 = %+d x_3 = %+d: bound %.5f, true smallest conductor %d (ratio %.5f); character values matched at %d of %d primes' % (
                'odd ' if kappa else 'even', x2, x3, np.exp(-t), true, np.exp(-t)/true, match, len(P.primes)), flush=True)
