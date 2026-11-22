"""Own segmented sieve (independent of the repo)."""
import numpy as np


def small_primes(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0]


_PR = {}


def mask(L, W):
    """bool array, M[i] = (L + i) is prime, for L >= 2."""
    lim = int((L + W) ** 0.5) + 1
    if lim not in _PR:
        _PR[lim] = small_primes(lim)
    M = np.ones(W, dtype=bool)
    for p in _PR[lim]:
        p = int(p)
        st = max(p * p, ((L + p - 1) // p) * p)
        if st >= L + W:
            continue
        M[st - L::p] = False
    if L <= 1:
        M[:2 - L] = False
    return M


def primes(L, W, seg=10 ** 7):
    out = []
    for s in range(L, L + W, seg):
        w = min(seg, L + W - s)
        out.append(np.nonzero(mask(s, w))[0].astype(np.int64) + s)
    return np.concatenate(out)
