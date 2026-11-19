"""Gilbreath depth against the largest gap: primes and fakes (true primes up to 100, then numbers without prime factors
up to B, chance prod p/(p-1) / ln n), 100 runs per B.  For each model: depth median, share of fakes at or below the
primes' depth, largest-gap median, share with largest gap <= the primes', and among those the share at or below the
primes' depth.  Prediction (from an independent review, before this run): the depth of the primes is typical among
fakes whose largest gap is as small as theirs.   usage: gilb4.py X runs out.txt"""
import sys
import numpy as np
X, RUNS = int(sys.argv[1]), int(sys.argv[2]); out = open(sys.argv[3], 'w')
def pr(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')
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
dp, _ = depth(P); gp = int(np.diff(P).max())
pr('X = %d: primes depth %d, largest gap %d' % (X, dp, gp))
rng = np.random.default_rng(8); m = np.arange(X + 1, dtype=float)
for B in (2, 7, 13, 47, 100):
    W = P[P <= B]; ok = np.ones(X + 1, bool); ok[:101] = False
    for p in W: ok[::p] = False
    c = float(np.prod(W / (W - 1.0))); pb = np.zeros(X + 1); pb[ok] = np.minimum(1, c / np.log(m[ok]))
    D = []; G = []
    for r in range(RUNS):
        s = np.concatenate([P[P <= 100], np.nonzero(rng.random(X + 1) < pb)[0]]); D.append(depth(s)[0]); G.append(int(np.diff(s).max()))
    D = np.array(D); G = np.array(G); sm = G <= gp
    pr('B = %3d: depth median %5.1f, at or below the primes %.2f; largest gap median %5.1f, <= primes %3d of %d; among those at or below the primes %s' % (B, np.median(D), np.mean(D <= dp), np.median(G), sm.sum(), RUNS, '%.2f' % np.mean(D[sm] <= dp) if sm.sum() else '-'))
