import numpy as np
from common import *
d = real_zeros(); Ds = sorted(d)
ls = np.round(np.arange(0.25, 8.01, 0.125), 3)
S = []; Ls = []
for D in Ds:
    s, L, n = unfold_stats(d[D], abs(D), 100, 600, ls); S.append(s); Ls.append(L)
S = np.array(S); Ls = np.array(Ls)
I = np.array([[1.0 if abs(D) % p == 0 else 0.0 for p in PR] for D in Ds])
X = np.column_stack([np.ones(len(Ds)), I])
print(' l     fit2   pred2  resid2  se | chi2(6 primes) | resid sum over primes')
for j, l in enumerate(ls):
    rest = np.array([sum(dS(l, p, Ls[k]) for p in primes_of(D) if p > 13) for k, D in enumerate(Ds)])
    b, se = ols(X, S[:, j] - rest)
    P = np.array([np.mean([dS(l, p, Ls[k]) for k in range(len(Ds)) if I[k, i]]) for i, p in enumerate(PR)])
    r = b[1:] - P
    print(f'{l:5.3f} {b[1]:.4f} {P[0]:.4f} {r[0]:+.4f} {se[1]:.4f} | {np.sum((r/se[1:])**2):6.1f} | {r.sum():+.4f}')
