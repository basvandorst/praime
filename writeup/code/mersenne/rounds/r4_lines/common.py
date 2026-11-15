"""Shared data: the Mersenne exponents, the bound B, the Wagstaff null.

Null (Wagstaff): every prime p < B is an exponent independently with
probability P(p) = min(1, e^gamma * log2(a*p) / p), a = 6 if p = 1 mod 4,
a = 2 if p = 3 mod 4; p = 2 is always included (it is an exponent and the
formula does not apply to it).

Exact sampler: for primes with P < 1 put mu = -log(1 - P); draw a Poisson
number M ~ Poisson(sum mu) of indices from the categorical law mu / sum mu
and keep the distinct ones.  A prime is hit at least once with probability
1 - exp(-mu) = P, independently of the others (Poisson splitting), so the
kept set is exactly the independent-Bernoulli null.  Conditioning on the
set size is done by rejection.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EXPS = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203,
        2281, 3217, 4253, 4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497,
        86243, 110503, 132049, 216091, 756839, 859433, 1257787, 1398269,
        2976221, 3021377, 6972593, 13466917, 20996011, 24036583, 25964951,
        30402457, 32582657, 37156667, 42643801, 43112609, 57885161, 74207281,
        77232917, 82589933]
B = 83880077
PAGE = 500000
assert len(EXPS) == 51 and EXPS == sorted(EXPS)
EULER_GAMMA = 0.57721566490153286061


def primes_below(n):
    """Odd-only numpy sieve; returns all primes < n as int64."""
    m = (n - 1) // 2          # index i <-> 2i+1, i = 0..m-1 covers odd < n
    s = np.ones(m, dtype=bool)
    s[0] = False              # 1
    r = int(n ** 0.5) + 1
    for i in range(1, (r - 1) // 2 + 1):
        if s[i]:
            p = 2 * i + 1
            s[(p * p - 1) // 2::p] = False
    odd = 2 * np.nonzero(s)[0].astype(np.int64) + 1
    return np.concatenate([np.array([2], dtype=np.int64), odd])


def wagstaff_prob(p):
    p = np.asarray(p, dtype=np.int64)
    a = np.where(p % 4 == 1, 6.0, 2.0)
    pr = np.exp(EULER_GAMMA) * np.log2(a * p) / p
    pr = np.minimum(pr, 1.0)
    pr[p == 2] = 1.0
    return pr


def load_primes():
    f = os.path.join(HERE, "primes_prob.npz")
    if os.path.exists(f):
        z = np.load(f)
        return z["p"].astype(np.int64), z["pr"]
    p = primes_below(B)
    pr = wagstaff_prob(p)
    np.savez(f, p=p.astype(np.uint32), pr=pr)
    return p, pr
