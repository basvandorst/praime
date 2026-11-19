"""Correlation of prime gaps at lags 1 to 6, and the variance of the turtle's sums g1 - g3 + g5 - ..., against the
odd Cramer model (random odd numbers, chance 2 / ln n).   usage: gapcorr.py X"""
import sys
import numpy as np
X = int(sys.argv[1])
sv = np.ones(X + 1, bool); sv[:2] = False
for i in range(2, int(X ** .5) + 1):
    if sv[i]: sv[i * i::i] = False
def stats(pr):
    P = np.nonzero(pr)[0]; P = P[P > 1000]
    g = np.diff(P).astype(float); g = g / np.log(P[1:])           # gaps scaled by the local mean gap
    g = g - g.mean()
    out = [float(np.mean(g[:-L] * g[L:]) / np.var(g)) for L in range(1, 7)]
    return out
rng = np.random.default_rng(3)
n = np.arange(X + 1, dtype=float); p2 = np.zeros(X + 1); p2[3::2] = 2 / np.log(n[3::2])
print('lag correlations of scaled gaps, lags 1..6')
print('primes      ', ' '.join('%+.4f' % v for v in stats(sv)))
v = np.array([stats(rng.random(X + 1) < p2) for _ in range(10)])
print('odd Cramer  ', ' '.join('%+.4f' % v for v in v.mean(0)), ' (spread', ' '.join('%.4f' % s for s in v.std(0)) + ')')
