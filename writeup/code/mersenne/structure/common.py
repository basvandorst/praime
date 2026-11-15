"""Shared data and segmented sieve for the Mersenne-exponent feature study."""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

EXPS = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281,
        3217, 4253, 4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243,
        110503, 132049, 216091, 756839, 859433, 1257787, 1398269, 2976221, 3021377,
        6972593, 13466917, 20996011, 24036583, 25964951, 30402457, 32582657,
        37156667, 42643801, 43112609, 57885161, 74207281, 77232917, 82589933]
assert len(EXPS) == 51

# GIMPS milestone page (fetched 2026-11-15): "All exponents below 83880077 have
# been tested and verified."  So the list is complete for p < B.
B = 83880077
SPLIT = 10**5          # discovery: p < SPLIT ; holdout: SPLIT <= p < B
EULER_GAMMA = 0.5772156649015329


def small_primes(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if s[i]:
            s[i*i::i] = False
    return np.flatnonzero(s)


def segmented_primes(N, seg=1 << 23):
    """All primes <= N as uint32, segmented over odd numbers (seg odd numbers per block)."""
    base = small_primes(int(N**0.5) + 1)
    base_odd = base[base > 2]
    out = [np.array([2], dtype=np.uint32)] if N >= 2 else []
    # block covers odd numbers lo, lo+2, ..., lo+2*(seg-1)
    lo = 3
    while lo <= N:
        hi = min(lo + 2 * seg, N + 1)          # exclusive
        n = (hi - lo + 1) // 2
        isp = np.ones(n, dtype=bool)
        for q in base_odd:
            q = int(q)
            if q * q >= hi:
                break
            start = max(q * q, ((lo + q - 1) // q) * q)
            if start % 2 == 0:
                start += q
            isp[(start - lo) // 2::q] = False
        vals = lo + 2 * np.flatnonzero(isp)
        vals = vals[vals <= N]
        out.append(vals.astype(np.uint32))
        lo = hi if hi % 2 == 1 else hi + 1
    return np.concatenate(out)


def load_primes():
    """Primes <= 2B+1 (cached)."""
    path = os.path.join(HERE, "primes_2B.npy")
    if os.path.exists(path):
        return np.load(path)
    pr = segmented_primes(2 * B + 1)
    np.save(path, pr)
    return pr


def null_prob(p, a_mode="mod4"):
    """Wagstaff/LPW heuristic P(2^p-1 prime) = e^gamma log2(a p)/p, capped at 1.
    a_mode='mod4': a=2 for p=3 mod 4 (and p=2), a=6 for p=1 mod 4.
    a_mode='a2'  : a=2 for all p."""
    p = np.asarray(p, dtype=np.float64)
    if a_mode == "mod4":
        a = np.where(np.mod(p, 4) == 1, 6.0, 2.0)
    else:
        a = np.full_like(p, 2.0)
    P = np.exp(EULER_GAMMA) * np.log2(a * p) / p
    return np.minimum(P, 1.0)
