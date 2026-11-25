"""Segmented odd-only sieve (numpy). Shared by tz_counts.py and the checks.

odd_flags(s, n, sp): boolean array f of length n with f[i] = (s + 2i is prime), s odd.
"""
import numpy as np

_SMALL = (3, 5, 7, 11, 13)
_PERIOD = 3 * 5 * 7 * 11 * 13          # period of the pattern in odd-index space


def base_primes(limit):
    """All primes <= limit (simple sieve, limit up to a few 1e6)."""
    a = np.ones(limit + 1, dtype=bool)
    a[:2] = False
    a[4::2] = False
    for p in range(3, int(limit ** 0.5) + 1, 2):
        if a[p]:
            a[p * p::2 * p] = False
    return np.flatnonzero(a)


def _pattern():
    # pattern[i] = False if 1 + 2i is divisible by a small prime (index space of odd numbers from 1)
    pat = np.ones(_PERIOD, dtype=bool)
    for p in _SMALL:
        # odd m = 1 + 2i divisible by p  <=> i = (p - 1)/2 mod p
        pat[(p - 1) // 2::p] = False
    return pat


_PAT = _pattern()


def odd_flags(s, n, sp):
    """Primality flags of the odd numbers s, s+2, ..., s+2(n-1). s must be odd, s >= 1.
    sp: array of odd primes covering sqrt(s + 2n)."""
    assert s % 2 == 1
    i0 = (s - 1) // 2                     # index of s in the odd-number space starting at 1
    off = i0 % _PERIOD
    reps = (off + n) // _PERIOD + 1
    f = np.tile(_PAT, reps)[off:off + n].copy()
    top = s + 2 * (n - 1)
    # the small primes themselves are prime
    for p in _SMALL:
        if s <= p <= top:
            f[(p - s) // 2] = True
    lim = int(top ** 0.5) + 1
    for p in sp:
        if p > lim:
            break
        if p <= 13:
            continue
        p = int(p)
        m = max(p * p, ((s + p - 1) // p) * p)
        if m % 2 == 0:
            m += p
        if m > top:
            continue
        f[(m - s) // 2::p] = False
    if s == 1:
        f[0] = False                      # 1 is not prime
    return f
