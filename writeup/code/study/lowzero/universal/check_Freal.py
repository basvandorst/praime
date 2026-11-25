"""Is F_p(X) = sum_k conj(lambda(p^k)) mu(p^k) X^k real for unitary Satake parameters (any degree)? and B_p = 1 + O(X^2)?"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from ulib import series_inv
rng = np.random.default_rng(1)
for d in (1, 2, 3, 4, 5):
    worst = 0; worstB = 0
    for trial in range(200):
        a = np.exp(2j * np.pi * rng.random(d))
        P = np.poly(a)[::1]          # coefficients highest first: prod (z - a) -> want prod(1 - a X)
        P = np.poly(a)               # [1, -e1, e2, ...] = prod(1 - a_i X) coefficients in X ascending order
        lam = series_inv(list(P), 60)
        X = 0.37
        F = sum(np.conj(lam[k]) * P[k] * X ** k for k in range(len(P)))
        worst = max(worst, abs(F.imag))
        # B_p = F * prod_{i,j} (1 - a_i conj(a_j) X)^-1 at small X: check 1 + O(X^2)
        for x in (1e-3, 1e-4):
            Fx = sum(np.conj(lam[k]) * P[k] * x ** k for k in range(len(P)))
            RS = np.prod([1 / (1 - ai * np.conj(aj) * x) for ai in a for aj in a])
            worstB = max(worstB, abs(Fx * RS - 1) / x ** 2)
    print('d=%d  max |Im F_p(0.37)| over 200 random unitary Satake sets: %.2e ;  max |B_p - 1|/X^2 at X=1e-3,1e-4: %.2f' % (d, worst, worstB))
