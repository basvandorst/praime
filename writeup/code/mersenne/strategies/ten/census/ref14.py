"""Independent reference for prime14.html: census of omega(n), n <= X (sympy primes, numpy), mean, variance, the three
predicted bar vectors for k = 0..11 (Erdos-Kac normal, Sathe-Selberg with G over the primes up to max(X, 10^6), the
normal with mean sum floor(X/p)/X and variance sum (1/p)(1-1/p) - sum_{p != q, pq > X} 1/(pq)), Poisson (Landau),
and the total variation distances over k = 0..11.   usage: ref14.py X out.json"""
import sys, json, math
import numpy as np
from sympy import primerange
from scipy.stats import norm
from scipy.special import gammaln
X = int(sys.argv[1]); K = 12
PG = np.array(list(primerange(2, max(X, 10 ** 6) + 1)), dtype=float); P = PG[PG <= X]
om = np.zeros(X + 1, np.int64)
for p in P.astype(np.int64): om[p::p] += 1
N = np.bincount(om[1:], minlength=K)[:K].astype(float)
L = math.log(math.log(X)); lx = math.log(X)
mean = float(np.sum(np.floor(X / P)) / X); v1 = float(np.sum((1 / P) * (1 - 1 / P)))
ex = 0.0
for p in P:
    ex += (1 / p) * float(np.sum(1 / P[P > X / p]))
    if p * p > X: ex -= 1 / (p * p)
v = v1 - ex
G = lambda z: math.exp(float(np.sum(np.log1p(z / (PG - 1)) + z * np.log1p(-1 / PG))) - gammaln(z + 1))
ek = [X * (norm.cdf((k + .5 - L) / math.sqrt(L)) - norm.cdf((k - .5 - L) / math.sqrt(L))) for k in range(K)]
tc = [X * (norm.cdf((k + .5 - mean) / math.sqrt(v)) - norm.cdf((k - .5 - mean) / math.sqrt(v))) for k in range(K)]
po = [0.0] + [X / lx * L ** (k - 1) / math.factorial(k - 1) for k in range(1, K)]
ss = [0.0] + [po[k] * G((k - 1) / L) for k in range(1, K)]
tv = lambda pv: float(np.sum(np.abs(N - np.array(pv))) / (2 * X))
mm = float(om[1:].mean()); vv = float(om[1:].var())
out = dict(N=N.astype(int).tolist(), mean=round(mm, 6), var=round(vv, 6), exmean=round(mean, 6), excl=round(ex, 6), tvar=round(v, 6), tv=[round(tv(ek), 6), round(tv(ss), 6), round(tv(tc), 6), round(tv(po), 6)])
json.dump(out, open(sys.argv[2], 'w')); print(out)
