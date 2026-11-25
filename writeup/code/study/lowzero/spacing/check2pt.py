"""Task 1(c), two points: does the wave closure reproduce the exact (Bogomolny-Keating, bklib) change of the number
variance when a prime is added to the Euler product?  closure: Sigma2_with(l) = E[Sigma2_base(l + D)] + Var(D) with the
base = the L-function with p missing (exact BK Sigma2), D the increment of p's wave over the window."""
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/lib')
import numpy as np
from scipy.interpolate import CubicSpline
import loadlibs, closure
bk = loadlibs.load_bklib(); tp = loadlibs.load_tplib(); G = tp.G
LG = bk.LGRID

def spline_even(vals):
    x = np.concatenate([[0.0], LG]); y = np.concatenate([[0.0], vals])
    cs = CubicSpline(np.concatenate([-x[:0:-1], x]), np.concatenate([y[:0:-1], y]))
    return lambda l: cs(np.abs(np.asarray(l, float)))

ells = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0]
idx = [int(round(l / bk.W)) - 1 for l in ells]
gue = bk.sigma2_gue()
for L in (9.4, 10.7):
    print('\nL = %.1f   rows: exact BK change (S -> S minus p, i.e. p made good) | closure, exact base | closure, sine base | Berry diagonal only' % L)
    print('  l       ' + ' '.join('%8.2f' % l for l in ells))
    for p, rest in ((2, ()), (3, ()), (5, ()), (7, ()), (2, (3,)), (3, (2, 5)), (5, (2,))):
        Sb = tuple(sorted(rest + (p,)))
        base = bk.sigma2_full(G, L, Sb); withp = bk.sigma2_full(G, L, tuple(rest))
        ex = (withp - base)[idx]
        w = closure.Wave(p, L)
        c_exact_base = closure.delta_count_var(spline_even(base), w, ells)
        c_sine = closure.delta_count_var(spline_even(gue), w, ells)
        berry = -bk.dsig_diag_grid(p, L)[idx]
        print('  p=%d rest=%-6s exact  ' % (p, str(rest)) + ' '.join('%+8.4f' % v for v in ex))
        print('                  closure' + ' '.join('%+8.4f' % v for v in c_exact_base))
        print('                  c-sine ' + ' '.join('%+8.4f' % v for v in c_sine))
        print('                  berry  ' + ' '.join('%+8.4f' % v for v in berry))
        print('                  ratio closure/exact ' + ' '.join('%8.3f' % (a / b) for a, b in zip(c_exact_base, ex)))
