import numpy as np, math
def base_primes(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]
def primes_in(lo, hi, bp=None):
    """primes in [lo, hi), lo odd-friendly, lo > sqrt(hi)"""
    if bp is None: bp = base_primes(int(math.isqrt(hi)) + 1)
    if lo % 2 == 0: lo0 = lo + 1
    else: lo0 = lo
    n = (hi - lo0 + 1) // 2           # odd numbers lo0, lo0+2, ...
    s = np.ones(n, bool)
    for p in bp[1:]:
        p = int(p)
        if p * p >= hi: break
        st = max(p * p, ((lo0 + p - 1) // p) * p)
        if st % 2 == 0: st += p
        s[(st - lo0) // 2::p] = False
    return lo0 + 2 * np.nonzero(s)[0].astype(np.int64)
