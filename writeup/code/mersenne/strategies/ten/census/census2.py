"""A fourth prediction for the census of omega(n): a normal law with the mean and variance of omega computed from the
primes themselves, with the exclusion that two primes p != q with pq > X cannot both divide n <= X:
  mean = sum_{p <= X} floor(X/p) / X,
  var  = sum_{p <= X} (1/p)(1 - 1/p) - sum_{p != q, pq > X} 1/(pq)        (heuristic; the last sum tends to pi^2/6).
Prediction written before the run: this normal law beats Erdos-Kac (mean = variance = ln ln X) by far, because the
exclusion term is large (about 1.6) at every reachable X.   usage: census2.py X out.txt"""
import sys, math
import numpy as np
from scipy.stats import norm
X = int(sys.argv[1]); out = open(sys.argv[2], 'w')
def pr(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')
sv = np.ones(X + 1, bool); sv[:2] = False
for i in range(2, int(X ** .5) + 1):
    if sv[i]: sv[i * i::i] = False
P = np.nonzero(sv)[0].astype(float)
om = np.zeros(X + 1, np.int8)
for p in P.astype(np.int64): om[p::p] += 1
N = np.bincount(om[1:]); K = len(N)
mean = float(np.sum(np.floor(X / P)) / X)
inv = 1 / P; cs = np.cumsum(inv[::-1])[::-1]
idx = np.searchsorted(P, X / P, side='right')
excl = float(np.sum(inv * np.where(idx < len(P), cs[np.minimum(idx, len(P) - 1)], 0.0)))
excl -= float(np.sum(inv * inv * (P * P > X)))
var = float(np.sum(inv * (1 - inv))) - excl
pr('X = %d: measured mean %.4f, variance %.4f; predicted mean %.4f, variance %.4f (exclusion term %.4f; pi^2/6 = %.4f); ln ln X = %.4f' % (X, float(om[1:].mean()), float(om[1:].var()), mean, var, excl, math.pi ** 2 / 6, math.log(math.log(X))))
pred = np.array([X * (norm.cdf((k + .5 - mean) / math.sqrt(var)) - norm.cdf((k - .5 - mean) / math.sqrt(var))) for k in range(K)])
pr('counts k = 0..%d measured: %s' % (K - 1, ' '.join(str(int(v)) for v in N)))
pr('normal with these constants:  %s' % ' '.join('%.0f' % v for v in pred))
pr('total variation distance %.4f' % float(np.sum(np.abs(N - pred)) / (2 * X)))
