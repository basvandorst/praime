"""Independent recount: twin pairs (p, p+2), p in [L, L+W), per class p mod 30030 (reduce to 2310, 210 later).
Odd-only segmented sieve, segments of 1e8 numbers. usage: count_twins.py L W out.npy"""
import sys, math
import numpy as np


def base_primes(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0]


def odd_mask(lo, hi, bp):
    """mask over odd numbers lo, lo+2, ..., < hi (lo odd): True if prime."""
    n = (hi - lo + 1) // 2
    m = np.ones(n, bool)
    for p in bp[1:]:
        p = int(p)
        if p * p >= hi:
            break
        st = max(p * p, ((lo + p - 1) // p) * p)
        if st % 2 == 0:
            st += p
        m[(st - lo) // 2::p] = False
    if lo <= 1:
        m[0] = False
    return m


L, W, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
Q = 30030
bp = base_primes(int(math.isqrt(L + W + 2)) + 2)
cnt = np.zeros(Q, np.int64)
seg = 10 ** 8
lo = L if L % 2 else L + 1
end = L + W
while lo < end:
    hi = min(lo + seg, end)
    m = odd_mask(lo, hi + 2, bp)            # includes p+2 for the last p < hi
    nlo = (hi - lo + 1) // 2                 # odd numbers in [lo, hi)
    t = np.nonzero(m[:nlo] & m[1:nlo + 1])[0]
    p = lo + 2 * t.astype(np.int64)
    cnt += np.bincount(p % Q, minlength=Q)
    lo = lo + 2 * nlo
np.save(out, cnt)
print(L, W, int(cnt.sum()))
