# Exploration after round 165: the measure at the first primes is forced to be a point mass (a genuine Satake angle).
# For every angle on the grid at p = 2 (and optionally a window at p = 3) the convex problem is solved for the rest;
# reported: the best t* over the grid.  Two forms (L(f), Sym^2).  Usage: pointmass.py <eps> <X> <r,r,...> [J]
import sys, time, numpy as np
from maass import Problem
eps = int(sys.argv[1]); X = float(sys.argv[2]); J = int(sys.argv[4]) if len(sys.argv) > 4 else 61
for r in [float(v) for v in sys.argv[3].split(',')]:
    t0 = time.time(); P = Problem(r, eps, X, kmax=2, J=J)
    lam0, tfree, _ = P.solve()
    best = (-9, None)
    for j in range(J):
        mask = np.ones((len(P.primes), J)); mask[0] = 0; mask[0, j] = 1
        lam, t, _ = P.solve(mask=mask)
        if t > best[0]: best = (t, j, P.summary())
    print('eps %d X %g r %.4f: free t* = %+.5f; point mass at 2: best t* = %+.5f at theta_2 = %.3f (a_2 = %+.3f) (%.0fs); a_p %s' % (
        eps, X, r, tfree, best[0], P.th[best[1]], 2*np.cos(P.th[best[1]]), time.time() - t0, best[2]), flush=True)
