# Controls with fixed data: the gap form must return the first zero.  zeta (with its pole), Delta, the character mod 3.
import numpy as np
from gap import Gap, tmax_fixed
tau = {2: -24, 3: 252, 5: 4830, 7: -16744, 11: 534612, 13: -577738, 17: -6905934, 19: 10661420, 23: 18643272, 29: 128406630, 31: -52843168}
def fixed(X, shifts, nu, logN, coef, pole=0):
    """coef(p, m) = c(p^m)/(d log p).  One column per prime power pattern: the measure is a point."""
    g = Gap(X, shifts, nth=3, pole=pole)
    g.cosm = np.array([[coef(p, m)] for p, m in g.pp]); sig = np.ones((len(g.primes), 1))
    return g, g.assemble(np.array(nu, float), logN, sig)
for X in (4.0, 8.0, 16.0, 32.0):
    g, H = fixed(X, [0.0], [1.0], 0.0, lambda p, m: 1.0, pole=1)
    print('zeta   X %g N %d: largest gap T %.4f (first zero 14.1347)' % (X, g.N, tmax_fixed(g, H)), flush=True)
    g, H = fixed(X, [5.5, 6.5], [0.5, 0.5], 0.0, lambda p, m: float(np.cos(m*np.arccos(tau[p]/(2*p**5.5)))))
    print('Delta  X %g N %d: largest gap T %.4f (first zero 9.2224)' % (X, g.N, tmax_fixed(g, H)), flush=True)
    chi = lambda p: 0 if p == 3 else (1 if p % 3 == 1 else -1)
    g, H = fixed(X, [1.0], [1.0], np.log(3), lambda p, m: float(chi(p))**m)
    print('chi_-3 X %g N %d: largest gap T %.4f (first zero 8.0397)' % (X, g.N, tmax_fixed(g, H)), flush=True)
