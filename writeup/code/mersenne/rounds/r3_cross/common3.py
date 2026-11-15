"""Shared data for the cross-domain feature factory (r3_cross).

Ranges (protocol of the coordinator):
  discovery  : odd primes 3 <= p < 1e5      (choose threshold + direction here)
  validation : 1e5 <= p < 1e7               (rank features, family-wise calibration)
  final exam : 1e7 <= p < B                 (only for 5 pre-registered features)
Null weight w(p) = e^g log2(a p)/p, a = 6 (p=1 mod 4), 2 (p=3 mod 4); w = 1 for p <= 7.
The segmented sieve is copied from writeup/code/mersenne/structure/common.py (read-only source).
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1"); os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EXPS = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281,
        3217, 4253, 4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243,
        110503, 132049, 216091, 756839, 859433, 1257787, 1398269, 2976221, 3021377,
        6972593, 13466917, 20996011, 24036583, 25964951, 30402457, 32582657,
        37156667, 42643801, 43112609, 57885161, 74207281, 77232917, 82589933]
B = 83880077
D_HI = 10**5
V_HI = 10**7
EG = 0.5772156649015329


def small_primes(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if s[i]:
            s[i*i::i] = False
    return np.flatnonzero(s)


def segmented_primes(N, seg=1 << 22):
    base = small_primes(int(N**0.5) + 1)
    base_odd = base[base > 2]
    out = [np.array([2], dtype=np.uint32)]
    lo = 3
    while lo <= N:
        hi = min(lo + 2 * seg, N + 1)
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
        out.append(vals[vals <= N].astype(np.uint32))
        lo = hi if hi % 2 == 1 else hi + 1
    return np.concatenate(out)


def primes_upto(N, tag):
    path = os.path.join(HERE, f"primes_{tag}.npy")
    if os.path.exists(path):
        return np.load(path)
    pr = segmented_primes(N)
    np.save(path, pr)
    return pr


def wnull(p):
    p = np.asarray(p, dtype=np.float64)
    a = np.where(np.mod(p, 4) == 1, 6.0, 2.0)
    w = np.exp(EG) * np.log2(a * p) / p
    return np.where(p <= 7, 1.0, np.minimum(w, 1.0))


# ---------- Poisson-binomial tails ----------
def pb_tail_from_powersums(c0, S, L=60):
    """Exact (to float) Poisson-binomial pmf from c0 = sum ln(1-w), S[k-1] = sum r^k (r=w/(1-w)).
    log G(z) = c0 + sum_k (-1)^(k+1) S_k z^k / k.  Returns tail[o] = P(X >= o), o = 0..L."""
    K = len(S)
    a = np.zeros(L + 1)
    for k in range(1, K + 1):
        if k <= L:
            a[k] = (-1) ** (k + 1) * S[k - 1] / k
    b = np.zeros(L + 1); b[0] = 1.0
    for n in range(1, L + 1):
        kk = np.arange(1, n + 1)
        b[n] = np.dot(kk * a[kk], b[n - kk]) / n
    pmf = np.exp(c0) * b
    pmf = np.maximum(pmf, 0)
    tail = np.cumsum(pmf[::-1])[::-1]
    return np.minimum(tail, 1.0)
