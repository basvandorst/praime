"""Measured statistics from count arrays mod 30030 (tcsieve.py output summed over the files of a window)."""
import math
import numpy as np
import fourpoint as F

def fold(a, Q):
    return np.bincount(np.arange(a.shape[-1]) % Q, weights=a, minlength=Q)

def measure(P, WP, X, Q, h=2):
    """P, WP, X: arrays of length 30030 (primes, sum 1/ln p, pairs (n, n+h) by n mod 30030)."""
    adm = F.admissible(Q, h); C = F.constants(Q, h); kap = C['S'] * F.phi(Q) / C['A']
    Xb = fold(X.astype(float), Q)[adm]; WPQ = fold(WP, Q)
    m = Xb.mean()
    Zb = kap * (WPQ[adm] + WPQ[[(b + h) % Q for b in adm]]); Zb = Zb - Zb.mean() + m
    Rb = Xb - Zb
    u2Z = float(np.var(Zb, ddof=1) / m); cov = float(np.cov(Xb, Zb)[0, 1] / m)
    return {'pairs': int(Xb.sum()), 'mean': float(m), 'u2X': float(np.var(Xb, ddof=1) / m), 'u2Z': u2Z,
            'u2R': float(np.var(Rb, ddof=1) / m), 'cov': cov, 'beta': cov / u2Z, 'K': len(adm)}
