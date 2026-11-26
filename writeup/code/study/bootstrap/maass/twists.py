# Exploration after round 165: add the twists f x chi and Sym^2 f x chi by real primitive characters chi_D (known data)
# to the Maass problem.  Level 1: conductor of f x chi is D^2, of Sym^2 f x chi is |D|^3.  Usage: twists.py <eps> <X> <r,r,..> <D,D,..|none> [kmax]
import sys, time, numpy as np
from sympy import kronecker_symbol
from maass import Problem, types
def twisted_types(r, eps, Ds, kmax=2):
    T = types(r, eps, kmax)
    for D in Ds:
        a = 0 if D > 0 else 1; e = (eps + a) % 2; q = abs(D)
        f1 = lambda ph: ph
        f1.pmp = (lambda D: lambda p, m, th: 2*np.cos(m*th)*int(kronecker_symbol(D, p))**m)(D)
        T.append(('f%d' % D, [e + 1j*r, e - 1j*r], 1, f1, 2*np.log(q)))
        f2 = lambda ph: ph
        f2.pmp = (lambda D: lambda p, m, th: (1 + 2*np.cos(2*m*th))*int(kronecker_symbol(D, p))**m)(D)
        T.append(('S2%d' % D, [a, a + 2j*r, a - 2j*r], 2, f2, 3*np.log(q)))
    return T
if __name__ == '__main__':
    eps = int(sys.argv[1]); X = float(sys.argv[2]); Ds = [] if sys.argv[4] == 'none' else [int(v) for v in sys.argv[4].split(',')]
    for r in [float(v) for v in sys.argv[3].split(',')]:
        t0 = time.time(); P = Problem(r, eps, X, T=twisted_types(r, eps, Ds), margin=12.0, N=max(16, int(np.ceil(np.log(X)/2*(2*r + 12)/np.pi))))
        lam, t, rd = P.solve()
        act = [x[0] for x, v in zip(P.T, P.mineigs(P.w)) if v < t + 1e-4]
        print('eps %d X %g r %.4f twists %s: t* = %.5f (%d rounds, %.0fs) active %s; a_p %s' % (eps, X, r, Ds, t, rd, time.time() - t0, act, P.summary()), flush=True)
