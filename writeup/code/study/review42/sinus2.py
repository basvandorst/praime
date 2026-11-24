import numpy as np, sys
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
from readings import *
from comparators import gls_chi2
for test in ('B', 'A'):
    S, full, diag, X, primes, ks = load(test)
    c2d, r, V = gls_chi2(S, diag, X); Vi = np.linalg.inv(V); npr = len(primes)
    def fitmodel(cols):
        A = np.zeros((16 * npr, len(cols) * npr))
        for j in range(16):
            for q in range(npr):
                A[j * npr + q, len(cols) * q:len(cols) * (q + 1)] = [f(prim[j]) for f in cols]
        beta = np.linalg.solve(A.T @ Vi @ A, A.T @ Vi @ r); rr = r - A @ beta
        return rr @ Vi @ rr, beta
    c = lambda l: np.cos(2 * np.pi * l)
    m1, _ = fitmodel([lambda l: c(l) * l])
    m2, b2 = fitmodel([c, lambda l: c(l) * l])
    m3, _ = fitmodel([c, lambda l: c(l) * l, lambda l: 1.0, lambda l: l])
    m4, _ = fitmodel([c, lambda l: c(l) * l, lambda l: c(l) * l * l, lambda l: 1.0, lambda l: l, lambda l: l * l])
    c2f, _, _ = gls_chi2(S, full, X)
    print('%s: GLS chi2 on %d cells: BK (0 free) %.1f | diag %.1f | diag + b*l*cos(2 pi l) per prime (1/prime) %.1f | (a + b l) cos (2/prime) %.1f | + c + d l (4/prime) %.1f | quadratic both (6/prime) %.1f' % (test, r.size, c2f, c2d, m1, m2, m3, m4))
