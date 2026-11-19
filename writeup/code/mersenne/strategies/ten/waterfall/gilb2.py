"""Gilbreath for fake sequences that copy the true primes up to B and are random (avoiding the multiples of W) after.
For each: depth, and the row and position where the first column breaks.   usage: gilb2.py X runs"""
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
rng = np.random.default_rng(2); n = np.arange(X + 1, dtype=float)
print('primes: depth %d, holds %s' % depth(P))
for B in (10, 100, 1000):
    for W in ((2,), (2, 3), (2, 3, 5)):
        ok = np.ones(X + 1, bool); ok[:B + 1] = False
        for p in W: ok[::p] = False
        c = 1.0
        for p in W: c *= p / (p - 1)
        pr = np.zeros(X + 1); pr[ok] = np.minimum(1, c / np.log(n[ok]))
        ds = []; br = 0
        for r in range(RUNS):
            s = np.concatenate([P[P <= B], np.nonzero(rng.random(X + 1) < pr)[0]])
            d, h = depth(s); ds.append(d); br += (not h)
        print('true primes up to %4d, then fakes without %-10s: depth mean %.1f, broken in %d of %d' % (B, W, np.mean(ds), br, RUNS))
