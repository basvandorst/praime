"""Independent reference for prime25.html: V(N) by direct gcd counting for small N (sympy-free, math.gcd) and by the phi
sieve for N up to 10^6; mean, spread and share positive of D(N)/N for N in [10^4, 10^6].   usage: ref25.py out.json"""
import sys, json, math
import numpy as np
X = 10 ** 6
phi = np.arange(X + 1, dtype=np.int64)
for p in range(2, X + 1):
    if phi[p] == p: phi[p::p] -= phi[p::p] // p
V = 2 * np.cumsum(phi) - 1
for n in (1, 2, 10, 24, 100):
    assert V[n] == sum(1 for a in range(1, n + 1) for b in range(1, n + 1) if math.gcd(a, b) == 1), n
C = 6 / math.pi ** 2; N = np.arange(10 ** 4, X + 1, dtype=float); d = (V[10 ** 4:] - C * N * N) / N
out = dict(V24=int(V[24]), V1e6=int(V[X]), mean=round(float(d.mean()), 9), sd=round(float(d.std()), 9), pos=round(float(np.mean(d > 0)), 9))
json.dump(out, open(sys.argv[1], 'w')); print(out)
