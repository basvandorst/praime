"""Smooth prediction of the number of trains of consecutive primes with sum <= x (so of the mean of f(n) over n <= x):
a train of k primes starting near t has sum about k t + (k^2/2) ln t (the gaps are about ln t); the starts have density
1/ln t, so the count of trains is about sum_{k >= 1} li(T_k), T_k the root of k t + (k^2/2) ln t = x, li(T) = int_2^T dt/ln t
(plus 1 for the start 2).  Compared with the exact count, and with Moser's limit ln 2.  Prediction written before the
run: the smooth count follows the exact one to within 1 %, and both lie below x ln 2 at every size tested.
usage: trains_model.py"""
import math
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
def li(T): return quad(lambda t: 1 / math.log(t), 2, T, limit=200)[0] if T > 2 else 0.0
def smooth(x):
    tot = 0.0; k = 1
    while True:
        g = lambda t: k * t + k * k / 2 * math.log(t) - x
        if g(2) >= 0: break
        T = brentq(g, 2, x); tot += li(T); k += 1
    return tot
def exact(x):
    sv = np.ones(x + 1, bool); sv[:2] = False
    for i in range(2, int(x ** .5) + 1):
        if sv[i]: sv[i * i::i] = False
    P = np.nonzero(sv)[0].astype(np.int64); cs = np.concatenate([[0], np.cumsum(P)]); n = 0
    for L in range(1, len(P) + 1):
        s = cs[L:] - cs[:-L]; c = int(np.sum(s <= x))
        if c == 0: break
        n += c
    return n
print('       x   exact trains  mean f   smooth model  ratio   x ln 2')
for x in (10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7):
    e = exact(x); s = smooth(x)
    print('%8d %12d %8.5f %14.1f %7.4f  %.5f' % (x, e, e / x, s, s / e, math.log(2)))
