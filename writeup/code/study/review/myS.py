"""Independent singular series S(B) = prod_p (1 - nu_p(B)/p) / (1 - 1/p)^k (Hardy-Littlewood), own implementation.
Primes up to X = 1e7 by a plain numpy sieve, tail sum_{p > X} [ln(1-k/p) - k ln(1-1/p)] ~ -C(k,2) / (X ln X)."""
import math
from functools import lru_cache
import numpy as np

X = 10_000_000


def _primes(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0]


PR = _primes(X).astype(np.float64)
PRI = PR.astype(np.int64)


@lru_cache(maxsize=None)
def logG(k, P):
    """sum over primes p > P of ln(1-k/p) - k ln(1-1/p), with p > k required (P >= k)."""
    p = PR[PR > P]
    assert P >= k - 1
    s = float(np.sum(np.log1p(-k / p) - k * np.log1p(-1.0 / p)))
    s += -k * (k - 1) / 2 / (X * math.log(X))
    return s


def S(B, P=None):
    B = sorted(set(int(b) for b in B))
    k = len(B)
    span = B[-1] - B[0]
    if P is None:
        P = max(span, k, 2)
    out = 0.0
    for p in PRI[PRI <= P]:
        p = int(p)
        nu = len(set(b % p for b in B))
        if nu >= p:
            return 0.0
        out += math.log(1 - nu / p) - k * math.log(1 - 1 / p)
    return math.exp(out + logG(k, P))


if __name__ == '__main__':
    print('S(0,2)     = %.10f' % S([0, 2]))
    print('2*C2 lit   = %.10f  (C2 = 0.6601618158468696)' % (2 * 0.6601618158468696))
    print('S(0,2,6)   = %.10f' % S([0, 2, 6]))
    print('S(0,4,6)   = %.10f' % S([0, 4, 6]))
    print('S(0,2,6,8) = %.10f' % S([0, 2, 6, 8]))
    print('S(0,6)     = %.10f' % S([0, 6]))
    print('S(0,30)    = %.10f' % S([0, 30]))
    print('S(0,2,4)   = %.10f' % S([0, 2, 4]))
    # Hardy-Littlewood triple constant: S(0,2,6) = (9/2) prod_{p>=5} (1-3/p)/(1-1/p)^3 * ... literature 2.858248596
    print('lit S(0,2,6) = 2.8582485957 (OEIS A065418 * 9/2: 0.6351663546 * 4.5 = %.10f)' % (0.6351663546 * 4.5))
