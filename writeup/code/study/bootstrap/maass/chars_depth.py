# Exploration (not registered): how many prescribed values make a real character extremal.  For each fundamental
# discriminant D and k = 0, 1, 2, ...: prescribe x_p = chi_D(p) at the first k primes, one measure on [-1, 1] at the others;
# bound(k) = smallest admissible conductor; truth(k) = smallest conductor of a real primitive character of the same parity
# that agrees with chi_D at those k primes.  Usage: chars_depth.py <X> <Dmax> <kmax>
import sys, numpy as np
from sympy import kronecker_symbol, primerange
from sympy.ntheory import factorint
from maass import Problem
X = float(sys.argv[1]); Dmax = int(sys.argv[2]); kmax = int(sys.argv[3])
def fundamental(D):
    if D % 4 == 1: return all(e == 1 for e in factorint(abs(D)).values())
    if D % 4 == 0:
        m = D//4; return m % 4 in (2, 3) and all(e == 1 for e in factorint(abs(m)).values())
    return False
grid = np.unique(np.concatenate([np.linspace(-1, 1, 81), np.linspace(-1, -0.98, 21), np.linspace(0.98, 1, 21), np.linspace(-0.02, 0.02, 21)]))
ps = list(primerange(2, 40)); fd = {s: [D for D in (range(-1, -400, -1) if s < 0 else range(2, 400)) if fundamental(D)] for s in (-1, 1)}
probs = {}
for s in (-1, 1):
    cf = lambda ph: ph; cf.pm = lambda m, x: x**m
    probs[s] = Problem(0.0, 0, X, T=[('chi', [1 if s < 0 else 0], 1, cf)], th=grid, margin=14.0)
for D in sorted(fd[-1][:60] + fd[1][:60], key=abs):
    if abs(D) > Dmax: break
    s = -1 if D < 0 else 1; P = probs[s]; out = []
    for k in range(kmax + 1):
        truth = next(abs(E) for E in fd[s] if all(kronecker_symbol(E, p) == kronecker_symbol(D, p) for p in ps[:k]))
        mask = np.ones((len(P.primes), len(grid)))
        for i in range(k): mask[i] = (np.abs(grid - int(kronecker_symbol(D, ps[i]))) < 1e-12)
        try: lam, t, rd = P.solve(mask=mask, tol=1e-7, rounds=150); b = np.exp(-t)
        except RuntimeError: b = float('nan')
        out.append((b, truth))
        if truth == abs(D) and b > 0.9999*abs(D): break
    depth = next((k for k, (b, tr) in enumerate(out) if tr == abs(D) and b > 0.9999*abs(D)), None)
    print('D = %4d: %s  -> exact from k = %s' % (D, '  '.join('k=%d: %.3f/%d' % (k, b, tr) for k, (b, tr) in enumerate(out)), depth), flush=True)
