"""The doorman: for each odd prime p, n(p) = the least quadratic non-residue mod p and g(p) = the least primitive root.
Erdos (1961): the share of p with n(p) = p_k tends to 2^-k, so the mean of n(p) tends to sum_k p_k / 2^k.
Elliott and Murata (1997): conjectured mean of g(p).  Prediction written before the run: shares 1/2, 1/4, 1/8 within
1 %, the mean of n(p) near 3.67 but approached slowly from below (rare large n(p) still missing).  usage: door.py X out.txt"""
import sys, math
import numpy as np
from sympy import primerange
X = int(sys.argv[1]); out = open(sys.argv[2], 'w')
def pr(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')
Q = list(primerange(2, 400))
E = sum(q / 2 ** (k + 1) for k, q in enumerate(Q))
pr('Erdos constant sum p_k / 2^k = %.6f' % E)
spf = np.zeros(X + 1, np.int64)
for q in primerange(2, int(X ** .5) + 2):
    m = spf[q * q::q]; m[m == 0] = q; spf[q * q::q] = m
P = [p for p in range(3, X + 1) if spf[p] == 0]
def pf(n):
    s = set()
    while n > 1:
        q = int(spf[n]) or n; s.add(q)
        while n % q == 0: n //= q
    return s
cnt = {}; S = 0; G = 0; rec = (0, 0); grec = (0, 0)
for p in P:
    h = (p - 1) // 2
    for q in Q:
        if pow(q, h, p) == p - 1: break
    cnt[q] = cnt.get(q, 0) + 1; S += q
    if q > rec[0]: rec = (q, p)
    F = pf(p - 1); g = 2
    while any(pow(g, (p - 1) // r, p) == 1 for r in F): g += 1
    G += g
    if g > grec[0]: grec = (g, p)
n = len(P)
pr('X = %d, %d odd primes' % (X, n))
pr('least non-residue: ' + ', '.join('%d: %.5f (2^-%d = %.5f)' % (q, cnt[q] / n, k + 1, 2 ** -(k + 1)) for k, q in enumerate(Q[:7]) if q in cnt))
pr('mean n(p) %.5f (Erdos %.5f); largest n(p) %d at p = %d' % (S / n, E, rec[0], rec[1]))
pr('mean least primitive root g(p) %.5f; largest %d at p = %d' % (G / n, grec[0], grec[1]))
