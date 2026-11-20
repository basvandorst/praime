"""Fake needles averaged over realisations: lattice points (a, b), a > 0, b >= 0, gcd(a, b) = 1, a^2 + b^2 = n odd, n <= X,
each kept with chance c / ln n (c gives the expected count N of the prime ideals).  For K = 16 ... the mean and spread
over R realisations of variance / mean of the sector counts (64 sector positions, as rw.py), next to the primes and the
Bernoulli floor 1 - (mean chance, weighted like the counts).  usage: rw_fakes.py X R out.txt"""
import sys, math
import numpy as np
X, RR = int(sys.argv[1]), int(sys.argv[2]); out = open(sys.argv[3], 'w')
def pr(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')
sv = np.ones(X + 1, bool); sv[:2] = False
for i in range(2, int(X ** .5) + 1):
    if sv[i]: sv[i * i::i] = False
A = []; LN = []; LT = []
for a in range(1, int(X ** .5) + 1):
    b = np.arange(0, int(math.isqrt(X - a * a)) + 1); n = a * a + b * b; t = np.arctan2(b, a)
    A.append(t[sv[n]]); ok = (np.gcd(a, b) == 1) & (n % 2 == 1) & (n >= 3); LN.append(n[ok]); LT.append(t[ok])
inert = [p for p in range(3, int(X ** .5) + 1, 4) if sv[p]]
A = np.concatenate(A + [np.zeros(len(inert))]); N = len(A); LN = np.concatenate(LN).astype(float); LT = np.concatenate(LT)
c = N / float(np.sum(1 / np.log(LN))); p = np.minimum(1, c / np.log(LN)); floor = 1 - float(np.sum(p * p) / np.sum(p))
def nv(t, K, offs=64):
    w = (math.pi / 2) / K; v = []
    for o in range(offs):
        cc = np.bincount(np.floor(((t + o * w / offs) % (math.pi / 2)) / w).astype(int), minlength=K)[:K]; v.append(cc.var())
    return float(np.mean(v)) / (len(t) / K)
rng = np.random.default_rng(77); F = [LT[rng.random(len(LN)) < p] for _ in range(RR)]
pr('X = %d, N = %d; Bernoulli floor of the fakes 1 - sum p^2 / sum p = %.3f; %d realisations' % (X, N, floor, RR))
pr('     K  logK/logN  primes   fakes mean (spread)   RW')
K = 16
while K <= N // 4:
    vs = [nv(f, K) for f in F]; b = math.log(K) / math.log(N)
    pr('%6d %9.3f %7.3f %10.3f (%.3f) %7.3f' % (K, b, nv(A, K), np.mean(vs), np.std(vs), min(1, 2 * b))); K *= 4
