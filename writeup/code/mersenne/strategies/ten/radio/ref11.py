"""Independent reference for prime11.html (numpy FFT, exact sums from residue counts with complex numbers): for the
primes up to X, the station strengths for q = 1..12 (mean over a coprime to q of |F(a/q)| / N, F(theta) = sum e(p theta)),
and the noise levels: mean of |F_k|^2 / N over k in [0, L/2] farther than w = max(3, round(20 L / 2^21)) bins from
round(L a / q) for every a/q with q <= Q (Q = 12, 50), L = the least power of 2 >= 2 (X + 1).  usage: ref11.py X out.json"""
import sys, json, cmath, math
import numpy as np
from math import gcd
from sympy import primerange
X = int(sys.argv[1]); P = list(primerange(2, X + 1)); N = len(P)
peaks = []
for q in range(1, 13):
    cnt = [0] * q
    for p in P: cnt[p % q] += 1
    vals = [abs(sum(cnt[r] * cmath.exp(2j * math.pi * a * r / q) for r in range(q))) for a in range(q) if gcd(a, q) == 1]
    peaks.append(round(sum(vals) / len(vals) / N, 6))
L = 1
while L < 2 * (X + 1): L *= 2
ind = np.zeros(L); ind[P] = 1; A = np.abs(np.fft.rfft(ind)) ** 2 / N; H = L // 2 + 1
fl = {}
for Q in (12, 50):
    m = np.ones(H, bool); w = max(3, round(20 * L / 2097152))
    for q in range(1, Q + 1):
        for a in range(0, q // 2 + 1):
            if gcd(a, q) != 1: continue
            c = round(L * a / q); m[max(0, c - w):min(H - 1, c + w) + 1] = False
    fl[Q] = round(float(A[m].mean()), 6)
out = dict(N=N, peaks=peaks, floor12=fl[12], floor50=fl[50]); json.dump(out, open(sys.argv[2], 'w')); print(out)
