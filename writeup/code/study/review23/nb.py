"""Fairer Poisson comparator for round 83: Cramer model on odd numbers (each odd number prime w.p. r = 2/lam):
D_k/2 is negative binomial, P(D = 2m) = C(m-1, k-1) r^k (1-r)^(m-k), m >= k."""
import json, math
import numpy as np
from scipy import stats
exec(open('r83check.py').read().split('# class chain')[0])
r = 2 / lam
for k in (2, 3, 4):
    h = np.array(c['res'][str(k)]['hist_even_D'])[1:]; N = h.sum(); cnt = band(h.astype(float))
    m = Ds // 2; nb = np.where(m >= k, stats.nbinom.pmf(m - k, k, r), 0.0)
    print(k, 'NB chi2', chi2(cnt, N * band(nb / nb.sum())))
