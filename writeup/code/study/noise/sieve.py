"""Segmented odd-only sieve (numpy), no repo imports."""
import math
import numpy as np

_BASE = {}


def base_primes(n):
    if n in _BASE:
        return _BASE[n]
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, math.isqrt(n) + 1):
        if s[i]:
            s[i * i::i] = False
    out = np.nonzero(s)[0].astype(np.int64)
    _BASE[n] = out
    return out


def primes_in(lo, hi, seg=10_000_000):
    """All primes p with lo <= p < hi (lo > 2 assumed), as int64, yielded segment by segment."""
    bp = base_primes(math.isqrt(hi) + 1)[1:]                  # odd base primes
    a = lo | 1
    while a < hi:
        b = min(a + seg, hi)
        n = (b - a + 1) // 2                                   # odd numbers a, a+2, ..., < b
        m = np.ones(n, bool)
        first = np.maximum(bp * bp, ((a + bp - 1) // bp) * bp)
        first = np.where(first % 2 == 0, first + bp, first)
        idx = (first - a) // 2
        for p, i in zip(bp.tolist(), idx.tolist()):
            if i < n:
                m[i::p] = False
        P = a + 2 * np.nonzero(m)[0].astype(np.int64)
        yield P[P < hi]
        a = b if b % 2 == 1 else b + 1
