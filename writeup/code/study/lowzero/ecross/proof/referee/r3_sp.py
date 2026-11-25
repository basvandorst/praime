"""Referee check 3 (own code): the stationary-phase constant of Theorem 2.
For one n with t_n* = 2 pi sqrt(n/N) = T (T = 250, 1000, 4000) and a smooth Gaussian weight H of width sigma = T/8:
  exact   A = (1/2pi) n^{-1/2} int H(t) e^{i(2theta(t) - t ln n)} dt     (trapezoid; theta by mpmath)
  simple  B = (e^{3 pi i/4}/sqrt 2) N^{-1/4} n^{-1/4} e(-2 sqrt(n/N)) H(t_n*)
  second  B2 = B * [1 + (i t/4)(H''/H)(t_n*) - i/(48 t) - i/(6t)]  (H'' term, the constant-weight correction, and the
               -1/(6t) phase of Stirling), all evaluated at t_n* (first order in 1/t)
Expect |A/B - 1| ~ c/T and |A/B2 - 1| = O(1/T^2)."""
import os, math, cmath
import numpy as np, mpmath as mp
mp.mp.dps = 30
N = 373.0; PI = math.pi
for T in (250.0, 1000.0, 4000.0):
    n = round(N * T * T / (4 * PI * PI)); ts = 2 * PI * math.sqrt(n / N)
    sig = T / 8
    Hf = lambda t: np.exp(-0.5 * ((t - ts) / sig) ** 2)
    h = 0.05
    tg = np.arange(ts - 7 * sig, ts + 7 * sig, h)
    lnN, l2p, lnn = mp.log(N), mp.log(2 * mp.pi), mp.log(n)
    ph = np.array([float(mp.fmod(t * lnN - 2 * t * l2p + 2 * mp.im(mp.loggamma(1 + 1j * mp.mpf(t))) - t * lnn, 2 * mp.pi)) for t in tg])
    A = (1 / (2 * PI)) * n ** -0.5 * np.sum(Hf(tg) * np.exp(1j * ph)) * h
    B = (cmath.exp(3j * PI / 4) / math.sqrt(2)) * N ** -0.25 * n ** -0.25 * cmath.exp(-4j * PI * math.sqrt(n / N)) * Hf(ts)
    H2overH = ((ts - ts) ** 2 / sig ** 4 - 1 / sig ** 2)
    B2 = B * (1 + 1j * ts / 4 * H2overH - 1j / (48 * ts) - 1j / (6 * ts))
    print('T=%6.0f n=%d: |A/B-1| = %.3e (x T = %.3f), |A/B2-1| = %.3e (x T^2 = %.2f)' % (T, n, abs(A / B - 1), abs(A / B - 1) * T, abs(A / B2 - 1), abs(A / B2 - 1) * T * T))
