"""Euler's fountain: count n in [0, N) with f(n) = n^2 + n + 41 prime, against the Bateman-Horn (Hardy-Littlewood
conjecture F) prediction  C * sum_{n<N} 1/ln f(n)  with  C = prod_p (1 - w(p)/p) / (1 - 1/p),  w(p) = number of roots of
f mod p = 1 + (-163/p) (Kronecker symbol; w(163) = 1).  Per prime the factor is 1 - chi(p)/(p - 1) (p != 163), and
prod_p (1 - chi(p)/p)^-1 = L(1, chi_-163) = pi/sqrt(163) (class number 1), so
C = (sqrt(163)/pi) * prod_p (1 - chi(p)/(p-1)) / (1 - chi(p)/p)   [absolutely convergent], p = 163 handled separately.
Primality: sieve the values f(n), n < N, by all primes p < sqrt(f(N)) through the roots of f mod p.
Prediction written before the run: the count is within 1 % of the Bateman-Horn value at N = 10^6 and 10^7.
usage: euler.py N out.txt"""
import sys, math
import numpy as np
from sympy import primerange, sqrt_mod, isprime
N = int(sys.argv[1]); out = open(sys.argv[2], 'w')
def pr(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')
def chi(p):                       # Kronecker symbol (-163 / p) for a prime p
    if p == 2: return -1          # -163 = 5 mod 8
    if p == 163: return 0
    return 1 if pow((-163) % p, (p - 1) // 2, p) == 1 else -1
# constant
lim = 10 ** 7
prod = 1.0
for p in primerange(2, lim):
    c = chi(p)
    if p == 163: prod *= (1 - 1 / p) / (1 - 1 / p) / (1 - 0 / p)   # w(163) = 1: factor (1 - 1/p)/(1 - 1/p) = 1; L-factor 1
    else: prod *= (1 - c / (p - 1)) / (1 - c / p)
C = math.sqrt(163) / math.pi * prod
pr('C (Bateman-Horn constant of n^2 + n + 41) = %.7f  (absolutely convergent product to %d)' % (C, lim))
# sieve
top = int(math.isqrt(N * N + N + 41)) + 1
alive = np.ones(N, bool)
n = np.arange(N, dtype=np.int64)
for p in primerange(2, top + 1):
    c = chi(p)
    if c == -1: continue
    # roots of n^2 + n + 41 = 0 mod p: n = (-1 +- sqrt(-163)) / 2
    if p == 2: continue
    r = sqrt_mod(-163 % p, p, all_roots=True)
    inv2 = (p + 1) // 2
    for s in set(r):
        a = ((-1 + s) * inv2) % p
        alive[a::p] = False
        if a < N and a * a + a + 41 == p: alive[a] = True       # f(a) is p itself
cnt = {}; bh = {}
logs = 1.0 / np.log((n * n + n + 41).astype(np.float64))
cs_alive = np.cumsum(alive); cs_bh = np.cumsum(logs)
pr('f(n) prime for n < 40: %d of 40 (Euler)' % int(cs_alive[39]))
pr('      N      primes f(n), n < N    Bateman-Horn C * sum 1/ln f(n)    ratio')
for e in range(2, int(math.log10(N)) + 1):
    k = 10 ** e
    pr('%8.0e %16d %26.1f %14.5f' % (k, int(cs_alive[k - 1]), C * cs_bh[k - 1], cs_alive[k - 1] / (C * cs_bh[k - 1])))
# spot check with sympy
bad = [m for m in list(range(0, 2000)) + list(range(N - 300, N)) if bool(alive[m]) != isprime(m * m + m + 41)]
pr('spot check against sympy.isprime on 2,300 values of n: %d disagreements' % len(bad))
