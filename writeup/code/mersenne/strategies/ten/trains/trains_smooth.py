"""Trains of a smooth prime sequence: q_m = li^{-1}(m) (real numbers, li(T) = int_2^T dt/ln t + li(2)-offset ignored),
m = 1, 2, ..., and the exact count of trains (runs of consecutive q's) with sum <= x, against the count for the true
primes.  Prediction written before the run: within 0.5 % of the primes at every x from 10^4 to 10^7, so the slow,
falling approach of the mean to ln 2 is a property of the density 1/ln t alone.   usage: trains_smooth.py"""
import math
import numpy as np
from scipy.special import expi
def li(T): return expi(np.log(T)) - expi(math.log(2))
def li_inv(m):
    T = np.maximum(3.0, m * np.log(np.maximum(m, 2)))          # Newton on li(T) = m
    for _ in range(40): T = np.maximum(2.0000001, T - (li(T) - m) * np.log(T))
    return T
def count(seq, x):
    cs = np.concatenate([[0.0], np.cumsum(seq)]); n = 0
    for L in range(1, len(seq) + 1):
        s = cs[L:] - cs[:-L]; c = int(np.sum(s <= x))
        if c == 0: break
        n += c
    return n
print('       x   primes: trains  mean    smooth q_m: trains  mean    ratio   (ln 2 = %.5f)' % math.log(2))
for x in (10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7):
    sv = np.ones(x + 1, bool); sv[:2] = False
    for i in range(2, int(x ** .5) + 1):
        if sv[i]: sv[i * i::i] = False
    P = np.nonzero(sv)[0].astype(float)
    M = int(li(np.array([x * 1.0]))[0]) + 2; q = li_inv(np.arange(1, M + 1, dtype=float)); q = q[q <= x]
    a = count(P, x); b = count(q, x)
    print('%8d %16d %7.4f %19d %7.4f %8.4f' % (x, a, a / x, b, b / x, b / a))
