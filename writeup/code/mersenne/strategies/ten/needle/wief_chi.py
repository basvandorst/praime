"""Second look at the evenness of q(p)/p: chi-square for 10 and 20 bins at x = 10^6, 3*10^6, 10^7, and the same for p
itself replaced by random uniform numbers (300 samples of the same size) to give the spread.  usage: wief_chi.py out.txt"""
import sys, math, gmpy2
import numpy as np
from sympy import primerange
out = open(sys.argv[1], 'w')
def pr(s): print(s); out.write(s + '\n')
P = np.array(list(primerange(3, 10 ** 7 + 1)), dtype=np.int64)
U = np.array([((int(gmpy2.powmod(2, int(p) - 1, int(p) * int(p))) - 1) // int(p)) / int(p) for p in P])
rng = np.random.default_rng(5)
for x in (10 ** 6, 3 * 10 ** 6, 10 ** 7):
    u = U[P <= x]; n = len(u)
    for B in (10, 20):
        c = np.bincount(np.minimum((u * B).astype(int), B - 1), minlength=B); chi = float(((c - n / B) ** 2 / (n / B)).sum())
        sims = [float(((np.bincount((rng.random(n) * B).astype(int), minlength=B) - n / B) ** 2 / (n / B)).sum()) for _ in range(300)]
        pr('x = %d, %d bins: chi-square %.2f; share of 300 uniform samples below it: %.3f' % (x, B, chi, float(np.mean(np.array(sims) < chi))))
