# The same structured problem in degree 1: a real "character".  Unknown: for each prime a probability measure on
# x in [-1, 1] (x = chi(p)), weight <x^m> log p at p^m (one measure for all powers).  One form, Gamma_R(s + kappa),
# no pole.  Output: the smallest admissible conductor and the means of x.  Usage: chars.py <kappa> <X,X,...>
import sys, numpy as np
from maass import Problem
def run(kappa, X, J=81, N=None):
    cf = lambda ph: ph; cf.pm = lambda m, x: x**m
    P = Problem(0.0, 0, X, T=[('chi', [kappa], 1, cf)], th=np.linspace(-1, 1, J), N=N, margin=14.0)
    lam, t, rd = P.solve()
    return P, lam, t, rd
if __name__ == '__main__':
    kappa = int(sys.argv[1])
    for X in [float(v) for v in sys.argv[2].split(',')]:
        P, lam, t, rd = run(kappa, X)
        print('kappa %d X %g N %d: smallest conductor %.6f (t* = %.6f, %d rounds); mean x: %s' % (kappa, X, P.N, np.exp(-t), t, rd,
              ' '.join('%d:%+.3f' % (p, float(w @ P.th)) for p, w in zip(P.primes, P.w))), flush=True)
