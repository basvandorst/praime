"""Own helpers (independent of the repo): sieve, singular series."""
import math
import numpy as np


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.nonzero(s)[0]


_BASE = primes_upto(60000)


def seg_mask(lo, n):
    """boolean primality of lo..lo+n-1 (lo >= 2), base primes up to 60000 (valid to 3.6e9)."""
    assert lo + n < 60000 ** 2
    M = np.ones(n, dtype=bool)
    for p in _BASE:
        p = int(p)
        if p * p >= lo + n:
            break
        st = max(p * p, -(-lo // p) * p)
        M[st - lo::p] = False
    if lo < 2:
        M[:2 - lo] = False
    return M


def sing(B, P=2_000_000, prs=None):
    """direct singular series of the set B: product over p <= P plus a tail estimate."""
    B = sorted(set(B)); k = len(B)
    if prs is None:
        prs = primes_upto(P)
    s = 0.0
    for p in prs:
        p = int(p)
        nu = len({b % p for b in B})
        if nu >= p:
            return 0.0
        s += math.log1p(-nu / p) - k * math.log1p(-1 / p)
    # tail: log factor ~ -C(k,2)/p^2 for p > P
    s += -k * (k - 1) / 2 / (P * math.log(P))
    return math.exp(s)
