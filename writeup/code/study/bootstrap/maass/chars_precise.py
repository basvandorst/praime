# Real characters, tighter solve: the smallest conductor to more digits, with a finer grid near the atoms.  Usage: chars_precise.py
import numpy as np
from maass import Problem
for kappa, target in [(1, 3.0), (0, 5.0)]:
    for X in (8.0, 14.0, 24.0):
        cf = lambda ph: ph; cf.pm = lambda m, x: x**m
        grid = np.unique(np.concatenate([np.linspace(-1, 1, 81), np.linspace(-1, -0.98, 41), np.linspace(0.98, 1, 41), np.linspace(-0.02, 0.02, 41)]))
        P = Problem(0.0, 0, X, T=[('chi', [kappa], 1, cf)], th=grid, margin=14.0)
        try: lam, t, rd = P.solve(tol=1e-9, rounds=2000, ncut=6)
        except RuntimeError as e: print("kappa %d X %g: LP solver failed at this tolerance (%s)" % (kappa, X, str(e)[:40])); continue
        print('kappa %d X %g: smallest conductor %.10f (target %g; difference %.2e; LP gap %.1e; %d rounds); x_p: %s' % (kappa, X, np.exp(-t), target, np.exp(-t) - target, t - lam, rd,
              ' '.join('%d:%+.4f' % (p, float(w @ P.th)) for p, w in zip(P.primes, P.w))), flush=True)
