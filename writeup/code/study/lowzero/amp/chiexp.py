"""EXPLORATORY (post hoc, both samples known). Does the number variance depend on the character values chi_D(2),
chi_D(3) at good primes (BK and the ratios conjecture: no, only |chi(p)|^2 enters)? Regress S - full on the design plus
c2 = chi_D(2) (0 if 2 | D) and c3 = chi_D(3) (0 if 3 | D); chi2 of each coefficient vector over the 16 lengths with the
residual covariance, and its projection on the p = 2 off-diagonal pattern (in units of that pattern).
usage: chiexp.py -> chiexp_out.txt"""
import os, sys, math
import numpy as np
from sympy import jacobi_symbol
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'bk'))
from common import *

def kron(D, p):
    if D % p == 0: return 0
    if p == 2: return 1 if D % 8 in (1, 7) else -1
    return int(jacobi_symbol(D % p, p))

def run(nm, Ds, S, F, Dg):
    X = design(Ds); n = len(Ds)
    c2 = np.array([kron(D, 2) for D in Ds], float); c3 = np.array([kron(D, 3) for D in Ds], float)
    Xe = np.column_stack([X, c2, c3]); A = np.linalg.pinv(Xe.T @ Xe) @ Xe.T
    Y = S - F; b = A @ Y; res = Y - Xe @ b; C = res.T @ res / (n - Xe.shape[1]); Vi = np.linalg.inv(Xe.T @ Xe)
    Ax = np.linalg.pinv(X.T @ X) @ X.T; off2 = ((Ax @ F) - (Ax @ Dg))[1]
    for j, lab in ((-2, 'chi(2)'), (-1, 'chi(3)')):
        v = b[j]; Vv = Vi[j, j] * C; chi = float(v @ np.linalg.inv(Vv) @ v)
        proj = float(v @ off2 / (off2 @ off2)); sp = math.sqrt(off2 @ Vv @ off2) / (off2 @ off2)
        print('  %-5s %s: coefficient chi2 over 16 lengths %.1f (null mean 16); projection on the p=2 off pattern %+.3f +- %.3f; counts +1/-1: %d/%d' % (
            nm, lab, chi, proj, sp, int((Xe[:, j] == 1).sum()), int((Xe[:, j] == -1).sum())))

Ds, S, F, Dg = load_B(); run('B', Ds, S, F, Dg)
z = np.load(os.path.join(HERE, 'r120_100_600.npz')); run('R120', [int(d) for d in z['D']], z['S'], z['F'], z['Dg'])
