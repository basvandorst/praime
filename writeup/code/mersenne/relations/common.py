"""Shared data and null models for the r2_relations experiment.

Real sequence: the 52 known Mersenne exponents, treated as plain integers.

Null (a): binomial (= conditional Poisson) process in x = log2 a.
  a_1 = 2 fixed; 51 points uniform on [1, log2 a_52(real)] (a homogeneous Poisson
  process conditioned on 51 points, i.e. the empirical rate 51/26.02 per bit),
  a_n = nearest ODD integer to 2^x (choice: odd, like every real term but 2);
  collisions are resampled.
  Note: log a_n = log log M_p - log log 2, so this is "log log M_p Poisson".

Null (b): Wagstaff-type matched null.  Each prime p is taken independently with
  probability min(1, e^gamma * log2(2p) / p).  Primes below X0 = 10^7 are enumerated
  exactly; above X0 the accepted primes are drawn as an inhomogeneous Poisson
  process with intensity e^gamma*log2(2x)/(x ln x) (prime density x acceptance),
  each point rounded to the next prime.  The first 52 accepted primes form a
  sequence; sequences are kept only if |log2 a_52 - log2 a_52(real)| <= 1
  (growth matched by rejection).
"""
import numpy as np
import gmpy2

A = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217,
     4253, 4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243, 110503, 132049,
     216091, 756839, 859433, 1257787, 1398269, 2976221, 3021377, 6972593, 13466917,
     20996011, 24036583, 25964951, 30402457, 32582657, 37156667, 42643801, 43112609,
     57885161, 74207281, 77232917, 82589933, 136279841]
assert len(A) == 52 and A == sorted(set(A))
N = 52
EG = float(np.exp(np.euler_gamma))
L52 = float(np.log2(A[-1]))


def null_a(rng):
    while True:
        x = np.sort(rng.uniform(1.0, L52, size=N - 1))
        v = [2]
        for xi in x:
            y = 2.0 ** xi
            o = int(np.floor((y - 1) / 2.0 + 0.5)) * 2 + 1  # nearest odd integer
            if o < 3:
                o = 3
            v.append(o)
        if len(set(v)) == N and v == sorted(v):
            return v


_X0 = 10 ** 7
_PR = None
_PROB = None


def _primes_below(n):
    s = np.ones(n // 2, dtype=bool)  # s[i] <-> 2i+1
    s[0] = False
    r = int(n ** 0.5) + 1
    for i in range(1, r // 2 + 1):
        if s[i]:
            p = 2 * i + 1
            s[p * p // 2::p] = False
    return np.concatenate(([2], 2 * np.nonzero(s)[0] + 1)).astype(np.int64)


def _init_b():
    global _PR, _PROB
    if _PR is None:
        _PR = _primes_below(_X0)
        pf = _PR.astype(np.float64)
        _PROB = np.minimum(1.0, EG * np.log2(2 * pf) / pf)


def null_b(rng, tol=1.0):
    _init_b()
    while True:
        acc = _PR[rng.random(_PR.size) < _PROB]
        seq = [int(v) for v in acc[:N]]
        if len(seq) < N:
            # continue above X0 by thinning in u = ln x
            u = np.log(_X0)
            lam_max = EG / np.log(2) * (1 + np.log(2) / u)  # intensity in u, decreasing
            while len(seq) < N:
                u += rng.exponential(1.0 / lam_max)
                lam = EG / np.log(2) * (1 + np.log(2) / u)
                if rng.random() < lam / lam_max:
                    p = int(gmpy2.next_prime(int(np.exp(u))))
                    if p > seq[-1]:
                        seq.append(p)
        if abs(np.log2(seq[-1]) - L52) <= tol:
            return seq


def make_nulls(n, seed=12345):
    rng = np.random.default_rng(seed)
    na = [null_a(rng) for _ in range(n)]
    nb = [null_b(rng) for _ in range(n)]
    return na, nb


def emp_p(real, null, side):
    """Empirical p-value, (1 + #as-extreme)/(n+1).  side: 'hi', 'lo' or 'two'."""
    null = np.asarray(null, dtype=float)
    n = null.size
    hi = (1 + np.sum(null >= real)) / (n + 1)
    lo = (1 + np.sum(null <= real)) / (n + 1)
    if side == 'hi':
        return hi
    if side == 'lo':
        return lo
    return min(1.0, 2 * min(hi, lo))
