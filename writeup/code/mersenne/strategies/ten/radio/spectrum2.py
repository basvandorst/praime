"""Spectrum of the primes up to X (as spectrum.py), noise floor measured away from the peaks a/q for all q <= Q,
for growing Q (w = 20 bins around a/q for q <= 50, 2 bins beyond).  Idea after spectrum.py: with Q = 50 the floor of
the primes equals that of fakes avoiding the primes up to 13; but the singular series predicts the peaks for every q,
and the sieve of the primes makes peaks up to q near sqrt(X), so the primes' 'noise' at Q = 50 is partly peaks.
Prediction written before the run: as Q grows the floor of the primes drops below that of the fakes with B = 13.
usage: spectrum2.py X runs out.txt"""
import sys
import numpy as np
from math import gcd
X, RUNS = int(sys.argv[1]), int(sys.argv[2]); out = open(sys.argv[3], 'w')
def pr(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')
sv = np.ones(X + 1, bool); sv[:2] = False
for i in range(2, int(X ** .5) + 1):
    if sv[i]: sv[i * i::i] = False
P = np.nonzero(sv)[0]; L = 1 << (2 * X).bit_length(); H = L // 2 + 1
QS = [q for q in (12, 50, 150, 400) if q * q < 0.6 * L / 5]
masks = {}
mask = np.ones(H, bool); done = 0
for Q in QS:
    for q in range(done + 1, Q + 1):
        ww = 20 if q <= 50 else 2
        for a in range(0, q // 2 + 1):
            if gcd(a, q) != 1: continue
            c = int(round(L * a / q)); mask[max(0, c - ww):min(H, c + ww + 1)] = False
    done = Q; masks[Q] = mask.copy()
pr('X = %d, L = %d; noise bins left for Q = %s: %s' % (X, L, QS, [int(masks[Q].sum()) for Q in QS]))
def floors(ind):
    A = np.abs(np.fft.rfft(ind.astype(float), L)) ** 2; N = ind.sum()
    return [float(A[masks[Q]].mean() / N) for Q in QS]
fp = floors(sv); pr('primes                       : floor ' + '  '.join('Q=%d %.4f' % (Q, v) for Q, v in zip(QS, fp)))
rng = np.random.default_rng(6); m = np.arange(X + 1, dtype=float)
for B in (2, 13, 100):
    W = P[P <= B]; ok = np.ones(X + 1, bool); cut = max(B, 100); ok[:cut + 1] = False
    for p in W: ok[::p] = False
    c = float(np.prod(W / (W - 1.0))); pb = np.zeros(X + 1); pb[ok] = np.minimum(1, c / np.log(m[ok])); pb[P[P <= cut]] = 1
    fl = np.array([floors(rng.random(X + 1) < pb) for s in range(RUNS)]).mean(0)
    pr('fakes avoiding primes up to %3d: floor ' % B + '  '.join('Q=%d %.4f' % (Q, v) for Q, v in zip(QS, fl)) + '  | primes / fakes ' + ' '.join('%.3f' % (a / b) for a, b in zip(fp, fl)))
