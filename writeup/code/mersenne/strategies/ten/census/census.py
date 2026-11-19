"""Census of omega(n), the number of different primes dividing n, for 1 <= n <= X, against three predictions:
(a) Erdos-Kac: normal with mean and variance ln ln X (counts from the normal on [k - 1/2, k + 1/2]);
(b) Landau / Poisson: N_k ~ (X / ln X) (ln ln X)^(k-1) / (k-1)!  for k >= 1;
(c) Sathe-Selberg: N_k ~ (X / ln X) (ln ln X)^(k-1) / (k-1)! G((k-1)/ln ln X),
    G(z) = prod_p (1 + z/(p-1)) (1 - 1/p)^z / Gamma(z + 1)   (product over p <= 10^7).
Measures: total variation distance sum |N_k - pred_k| / (2X), and the ratio pred/N for each k.
Prediction written before the run: Sathe-Selberg is far closer than Poisson, and Poisson closer than Erdos-Kac.
usage: census.py X out.txt"""
import sys, math
import numpy as np
from scipy.special import gammaln
from scipy.stats import norm
X = int(sys.argv[1]); out = open(sys.argv[2], 'w')
def pr(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')
sv = np.ones(X + 1, bool); sv[:2] = False
for i in range(2, int(X ** .5) + 1):
    if sv[i]: sv[i * i::i] = False
P = np.nonzero(sv)[0]
om = np.zeros(X + 1, np.int8)
for p in P: om[p::p] += 1
N = np.bincount(om[1:]); K = len(N)
L = math.log(math.log(X)); lx = math.log(X)
PP = P if X >= 10 ** 7 else None
if PP is None:
    M = 10 ** 7; s2 = np.ones(M + 1, bool); s2[:2] = False
    for i in range(2, int(M ** .5) + 1):
        if s2[i]: s2[i * i::i] = False
    PP = np.nonzero(s2)[0]
lp = np.log1p(-1.0 / PP)
def G(z): return math.exp(float(np.sum(np.log1p(z / (PP - 1.0)) + z * lp)) - gammaln(z + 1))
ek = np.array([X * (norm.cdf((k + .5 - L) / math.sqrt(L)) - norm.cdf((k - .5 - L) / math.sqrt(L))) for k in range(K)])
po = np.array([0.0] + [X / lx * L ** (k - 1) / math.factorial(k - 1) for k in range(1, K)])
ss = np.array([0.0] + [X / lx * L ** (k - 1) / math.factorial(k - 1) * G((k - 1) / L) for k in range(1, K)])
pr('X = %d, ln ln X = %.4f; measured mean of omega %.4f, variance %.4f' % (X, L, float(om[1:].mean()), float(om[1:].var())))
pr(' k   measured   Erdos-Kac    Poisson  Sathe-Selberg   (ratios prediction / measured)')
for k in range(K):
    pr('%2d %10d %11.1f %10.1f %14.1f   %6.3f %6.3f %6.3f' % (k, N[k], ek[k], po[k], ss[k], ek[k] / N[k], (po[k] / N[k]) if k else float('nan'), (ss[k] / N[k]) if k else float('nan')))
for name, pv in (('Erdos-Kac', ek), ('Poisson (Landau)', po), ('Sathe-Selberg', ss)):
    pr('%-17s total variation distance %.4f' % (name, float(np.sum(np.abs(N - pv)) / (2 * X))))
pr('smallest n with omega = k: ' + ', '.join('%d: %d' % (k, int(np.argmax(om == k))) for k in range(1, K)))
