"""Independent reference for prime16.html (primes only): for each width h and every n in [N/2, N], the variance / mean
of pi(n + h) - pi(n) - (li(n + h) - li(n)) with li(m) = sum_{2 <= t <= m} 1/ln t, the variance of psi(n + h) - psi(n) - h,
the Montgomery-Soundararajan value h (ln(n/h) + B) averaged, Cramer's h ln n averaged.   usage: ref16.py N out.json"""
import sys, json, math
import numpy as np
from sympy import primerange
N = int(sys.argv[1]); HS = [10, 30, 100, 300, 1000, 3000, 10000]; M = N + 10001
P = list(primerange(2, M + 1)); ind = np.zeros(M + 1); ind[P] = 1
lam = np.zeros(M + 1)
for p in P:
    q = p
    while q <= M: lam[q] += math.log(p); q *= p
psi = np.cumsum(lam); cs = np.cumsum(ind)
li = np.zeros(M + 1); li[2:] = np.cumsum(1 / np.log(np.arange(2, M + 1, dtype=float)))
B = 1 - 0.5772156649015329 - math.log(2 * math.pi)
n = np.arange(N // 2, N + 1); out = {}
for h in HS:
    d = (cs[n + h] - cs[n]) - (li[n + h] - li[n]); mu = float((li[n + h] - li[n]).mean())
    w = psi[n + h] - psi[n] - h
    out[str(h)] = dict(ratio=round(float(d.var()) / mu, 6), wv=round(float(w.var()), 4), ms=round(float(np.mean(h * (np.log(n / h) + B))), 4), cr=round(float(np.mean(h * np.log(n))), 4))
json.dump(out, open(sys.argv[2], 'w')); print(out)
