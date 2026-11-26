# Up to which support is the empty system feasible (zeta datum)?  Bisection on the smallest eigenvalue of the
# pole-plus-archimedean form; Fourier basis, float64.  Second: the true primes below p_(n+1) without p_(n+1).
import numpy as np
from sympy import primerange, nextprime
from conductor import FormG
for N in (20, 40, 60):
    lo, hi = 2.0, 2.2
    for _ in range(30):
        mid = (lo + hi)/2
        if FormG(mid, N, [0]).mineig([], []) > 0: lo = mid
        else: hi = mid
    print('N=%d: empty system positive up to X = %.5f (L = %.5f); lam_min at X = 2: %.4e' % (N, lo, np.log(lo), FormG(2.0, N, [0]).mineig([], [])), flush=True)
for pn in (2, 3, 5, 7, 11):
    q = int(nextprime(pn)); P = [float(p) for p in primerange(2, pn + 1)]
    out = []
    for fac in (1.0, 1.02, 1.05):
        X = q*fac; N = max(12, int(8 + 4*np.log(X)/np.log(2)))
        out.append('%.2e' % FormG(X, N, [0]).mineig(P, [1]*len(P)))
    print('primes up to %d only: lam_min at X = %d * (1.00, 1.02, 1.05): %s' % (pn, q, ' '.join(out)), flush=True)
