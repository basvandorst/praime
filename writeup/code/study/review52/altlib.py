"""Alternatives to the three-point prediction, matched through the registered design on the fresh data.
 closure : wave closure (explicit-formula wave of each bad prime + independence), all orders, per L-function
           sum over its bad primes at its zero-weighted L.
 first   : first order -Li2/2 + (ln p/L) l(2+l)   (= LDM)
 M-only  : 3pt with Delta V3,P replaced by the far-field model M (uses derivation's L^2*residual table: NOT own)
 noT3    : 3pt with Delta V3 = 0 (GUE T3 or T3 without bad-prime factors) (uses derivation's tables: NOT own)
Also the mean shift under the closure, and the 'remove the wave' map on the real values."""
import os, numpy as np, json, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from myfit import load, stat, design, ols, PR
from mpmath import polylog
from functools import lru_cache


def factor(n):
    f = []; d = 2
    while d * d <= n:
        if n % d == 0:
            f.append(d)
            while n % d == 0: n //= d
        d += 1
    if n > 1: f.append(n)
    return f


@lru_cache(maxsize=None)
def clo(p, Lr, n=1 << 13):
    L = Lr / 1000.
    K = max(4, int(60 / np.log2(p)) + 2)
    th = 2 * np.pi * np.arange(n) / n
    k = np.arange(1, K + 1)[:, None]
    c = p ** (-k / 2.0) * np.cos(k * th)
    A = (c / k).sum(0)
    eps = -(2 * np.log(p) / L) * c.sum(0)
    w = 1 + eps
    X = A + np.log(w)
    m = np.mean(w * X)
    return -(np.mean(w * X * X) - m * m), -m


def first(p, L):
    l = np.log(p / (p - 1.0))
    return -0.5 * float(polylog(2, 1.0 / p)) + np.log(p) / L * l * (2 + l)


if __name__ == "__main__":
  pass
