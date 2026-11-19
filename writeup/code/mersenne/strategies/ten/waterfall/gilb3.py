"""Depth of Gilbreath's triangle: the primes against 100 fake sequences per model (true primes up to 100 first, then
random numbers avoiding the multiples of W, chance prod p/(p-1) / ln n). Percentile of the primes.  usage: gilb3.py X runs"""
import sys
import numpy as np
X, RUNS = int(sys.argv[1]), int(sys.argv[2])
sv = np.ones(X + 1, bool); sv[:2] = False
for i in range(2, int(X ** .5) + 1):
    if sv[i]: sv[i * i::i] = False
P = np.nonzero(sv)[0]
def depth(seq):
    a = np.array(seq, dtype=np.int64); d = 0
    while True:
        if d > 0:
            if a[0] != 1: return d, False
            if len(a) <= 1 or np.all((a[1:] == 0) | (a[1:] == 2)): return d, True
        if len(a) <= 1: return d, True
        a = np.abs(np.diff(a)); d += 1
rng = np.random.default_rng(4); n = np.arange(X + 1, dtype=float)
dp, hp = depth(P)
print('X = %d: primes depth %d (holds %s)' % (X, dp, hp))
for W in ((2,), (2, 3), (2, 3, 5), (2, 3, 5, 7), (2, 3, 5, 7, 11, 13)):
    ok = np.ones(X + 1, bool); ok[:101] = False
    for p in W: ok[::p] = False
    c = 1.0
    for p in W: c *= p / (p - 1)
    pr = np.zeros(X + 1); pr[ok] = np.minimum(1, c / np.log(n[ok]))
    ds = []; br = 0
    for r in range(RUNS):
        d, h = depth(np.concatenate([P[P <= 100], np.nonzero(rng.random(X + 1) < pr)[0]])); ds.append(d); br += (not h)
    ds = np.array(ds)
    print('fakes without %-20s depth mean %6.1f, median %5.1f, 10%%..90%% %3d..%3d, share at or below the primes %.2f, broken %d of %d' % (W, ds.mean(), np.median(ds), np.percentile(ds, 10), np.percentile(ds, 90), np.mean(ds <= dp), br, RUNS))
