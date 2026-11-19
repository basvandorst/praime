"""Gilbreath's triangle: depth D (the first row whose entries after the leading 1 are all 0 or 2; from there on every
row starts with 1) and whether the first column stays 1, for the primes and for fake sequences.
usage: gilb.py X runs"""
import sys
import numpy as np
X, RUNS = int(sys.argv[1]), int(sys.argv[2])
def depth(seq):
    a = np.array(seq, dtype=np.int64); d = 0
    while True:
        if d > 0:
            if a[0] != 1: return d, False, int(a.max())
            if len(a) <= 1 or np.all((a[1:] == 0) | (a[1:] == 2)): return d, True, int(a.max())
        if len(a) <= 1: return d, True, 0
        a = np.abs(np.diff(a)); d += 1
sv = np.ones(X + 1, bool); sv[:2] = False
for i in range(2, int(X ** .5) + 1):
    if sv[i]: sv[i * i::i] = False
P = np.nonzero(sv)[0]
d, ok, mx = depth(P); print('primes up to %d (%d terms): depth %d, first column all 1: %s' % (X, len(P), d, ok))
rng = np.random.default_rng(1); n = np.arange(X + 1, dtype=float)
def fake(W):
    ok = np.ones(X + 1, bool); ok[:3] = False
    for p in W: ok[::p] = False
    c = 1.0
    for p in W: c *= p / (p - 1)
    pr = np.zeros(X + 1); pr[ok] = np.minimum(1, c / np.log(n[ok]))
    s = np.nonzero(rng.random(X + 1) < pr)[0]
    head = [p for p in (2, 3, 5, 7, 11, 13) if p in W]
    return np.concatenate([np.array(head + ([2] if not head else []), dtype=np.int64), s[s > (max(W) if W else 2)]])
for W in ((), (2,), (2, 3), (2, 3, 5), (2, 3, 5, 7)):
    ds = []; fails = 0
    for r in range(RUNS):
        d, ok, mx = depth(fake(W)); ds.append(d); fails += (not ok)
    print('fakes without %-12s: depth mean %.1f (min %d, max %d), first column broken in %d of %d' % (W or '()', np.mean(ds), min(ds), max(ds), fails, RUNS))
