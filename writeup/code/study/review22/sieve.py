"""Own segmented odd-only sieve. usage: sieve.py L W out.npy  -> primes in [L, L+W) as int64."""
import sys, math
import numpy as np


def base_primes(n):
    s = np.ones(n // 2 + 1, bool); s[0] = False           # index i <-> 2i+1
    for i in range(1, (math.isqrt(n) - 1) // 2 + 1):
        if s[i]:
            p = 2 * i + 1
            s[(p * p) // 2::p] = False
    return np.concatenate([[2], 2 * np.nonzero(s)[0] + 1]).astype(np.int64)


def primes_in(L, W, seg=20_000_000):
    bp = base_primes(math.isqrt(L + W) + 1)[1:]
    out = []
    lo = L if L % 2 else L + 1                               # odd numbers lo, lo+2, ...
    hi_all = L + W
    while lo < hi_all:
        hi = min(lo + seg, hi_all)
        n = (hi - lo + 1) // 2
        m = np.ones(n, bool)
        for p in bp.tolist():
            if p * p >= hi:
                break
            st = max(p * p, -(-lo // p) * p)
            if st % 2 == 0:
                st += p
            m[(st - lo) // 2::p] = False
        v = lo + 2 * np.nonzero(m)[0].astype(np.int64)
        out.append(v[v < hi])
        lo = lo + 2 * n
    P = np.concatenate(out)
    if L <= 2 < L + W:
        P = np.concatenate([[2], P[P > 2]])
    return P[P >= max(L, 2)]


if __name__ == '__main__':
    L, W = int(float(sys.argv[1])), int(float(sys.argv[2]))
    P = primes_in(L, W)
    np.save(sys.argv[3], P)
    print(len(P), P[:3], P[-3:])
