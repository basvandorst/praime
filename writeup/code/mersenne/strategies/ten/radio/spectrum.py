"""Spectrum of the primes up to X: |F(k/L)|^2 with F(theta) = sum_{p <= X} e(p theta), L = 2^m >= 2X.
Peaks sit at a/q (q square-free) with |F| about pi(X) / phi(q). Noise floor: the mean of |F|^2 / N over the bins
farther than w bins from every a/q with q <= Q.  Parseval: sum_k |F_k|^2 = L N, so mass in the peaks leaves less noise.
Prediction written before the run: the floor of the primes is lower than that of fakes with the same density, by
about the factor seen in the width of the Goldbach comet (0.645 against odd fakes), and the width ratio follows the
floor ratio for the sieved fakes too.   usage: spectrum.py X Q w runs out.txt"""
import sys
import numpy as np
from math import gcd
X, Q, w, RUNS = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]); out = open(sys.argv[5], 'w')
def pr(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')
sv = np.ones(X + 1, bool); sv[:2] = False
for i in range(2, int(X ** .5) + 1):
    if sv[i]: sv[i * i::i] = False
P = np.nonzero(sv)[0]
L = 1 << (2 * X).bit_length()
mask = np.ones(L // 2 + 1, bool)
for q in range(1, Q + 1):
    for a in range(0, q + 1):
        if gcd(a, q) != 1 or 2 * a > q: continue
        c = int(round(L * a / q)); mask[max(0, c - w):min(L // 2 + 1, c + w + 1)] = False
pr('X = %d, L = %d, bins kept as noise: %d of %d (Q = %d, w = %d)' % (X, L, mask.sum(), L // 2 + 1, Q, w))
def floor(ind):
    F = np.fft.rfft(ind.astype(float), L); A = np.abs(F) ** 2; N = ind.sum()
    return float(A[mask].mean() / N), float(A[~mask].sum() / (A.sum()))
fp, mp = floor(sv)
pr('primes: noise floor %.4f (mean |F|^2 / N on the noise bins); share of the half-spectrum mass near a/q, q <= Q: %.4f' % (fp, mp))
rng = np.random.default_rng(5); m = np.arange(X + 1, dtype=float)
for B in (2, 13, 31, 100):
    W = P[P <= B]; ok = np.ones(X + 1, bool); cut = max(B, 100); ok[:cut + 1] = False
    for p in W: ok[::p] = False
    c = float(np.prod(W / (W - 1.0))); pb = np.zeros(X + 1); pb[ok] = np.minimum(1, c / np.log(m[ok])); pb[P[P <= cut]] = 1
    fl = []; ms = []
    for s in range(RUNS):
        a, b = floor(rng.random(X + 1) < pb); fl.append(a); ms.append(b)
    pr('fakes avoiding primes up to %3d: noise floor %.4f (spread %.4f), mass near a/q %.4f; primes / fakes floor %.3f' % (B, np.mean(fl), np.std(fl), np.mean(ms), fp / np.mean(fl)))
