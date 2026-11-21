"""Theory predictions for the probes, each with its formula.

lam(t) = ln t.  All window averages are taken over the primes of the window, i.e. with weight 1/ln t, and the
integrals are done in u = ln t (dt = e^u du), which avoids the trapezoid bias near a 1/ln t end point.

Hardy-Littlewood singular series of a set H = {h_1..h_k}:
    S(H) = prod_p (1 - nu_p(H)/p) / (1 - 1/p)^k,   nu_p(H) = number of distinct residues of H mod p.
    (conjecture: Hardy-Littlewood 1923; S({0,2}) = 2 C2, C2 = 0.6601618158...)

Crude consecutive-gap law ("Gallagher-Brent", heuristic): the share of primes p ~ t followed by the gap h is
    S({0,h}) exp(-h/lam) / lam           (pair singular series times a Poisson chance of no prime in between)
and two equal gaps h, h in a row:  S({0,h,2h}) exp(-2h/lam) / lam^2.

Montgomery-Soundararajan (conjecture, 2004): counts of primes in [n, n+H) around n ~ x have
    Var / mean = 1 - (ln H - B)/ln x,   B = 1 - gamma - ln(2 pi).

Chebyshev bias (Rubinstein-Sarnak 1994, under GRH + LI): the mean of (pi(x;q,N) - pi(x;q,R)) ln x / sqrt x is
    (c(q,R) - c(q,N)) / phi(q),   c(q,a) = #{b mod q : b^2 = a};  = 1 for (4;3,1), (3;2,1), (8;a,1).

Bateman-Horn (conjecture, 1962) for f(n) = a n^2 + b n + c irreducible:
    #{n in I : f(n) prime} ~ C(f) sum_{n in I} 1/ln f(n),  C(f) = prod_p (1 - w(p)/p)/(1 - 1/p),
    w(p) = number of roots of f mod p.
"""
import math
import numpy as np
from common import small_primes

EULER = 0.5772156649015329
B_MS = 1 - EULER - math.log(2 * math.pi)
_PR = small_primes(200_000)


from functools import lru_cache


def singular(H):
    return _singular(tuple(sorted(set(int(h) for h in H))))


@lru_cache(maxsize=None)
def _singular(H):
    H = list(H)
    k = len(H)
    s = 1.0
    for p in _PR:
        p = int(p)
        nu = len(set(h % p for h in H))
        if nu == p:
            return 0.0
        s *= (1 - nu / p) / (1 - 1 / p) ** k
        if p > 4 * max(H + [2]) and nu == k and p > 50000:
            # tail: (1-k/p)/(1-1/p)^k = 1 - C(k,2)/p^2 + ..., summed over p > P: ~ C(k,2)/(P ln P)
            P = p
            s *= math.exp(-k * (k - 1) / 2 / (P * math.log(P)))
            break
    return s


def _window_avg(L, W, fn, n=4000):
    """sum over primes in [L, L+W) of fn(lam) / number of primes, model density 1/lam, integral in u = ln t."""
    u = np.linspace(math.log(L), math.log(L + W), n + 1)
    w = np.exp(u) / u                    # dt/ln t = e^u du / u
    f = fn(u)
    num = np.trapezoid(w * f, u)
    den = np.trapezoid(w, u)
    return float(num / den)


def gap_share(h, L, W):
    S = singular([0, h])
    return _window_avg(L, W, lambda lam: S * np.exp(-h / lam) / lam)


def eqgap_share(L, W, hmax=200):
    tot = 0.0
    for h in range(2, hmax + 1, 2):
        S = singular([0, h, 2 * h])
        if S:
            tot += _window_avg(L, W, lambda lam, S=S, h=h: S * np.exp(-2 * h / lam) / lam ** 2)
    return tot


def same_class_share(q, L, W, hmax=600):
    tot, norm = 0.0, 0.0
    for h in range(2, hmax + 1, 2):
        v = gap_share(h, L, W)
        norm += v
        if h % q == 0:
            tot += v
    return tot / norm


def ms_ratio(H, L, W):
    return _window_avg(L, W, lambda lam: 1 - (math.log(H) - B_MS) / lam)


def bh_constant(a, b, c, pmax=2_000_000):
    """C(f) for f = a n^2 + b n + c; w(p) by brute force for p <= 2000, else 1 + (D/p) for p not dividing aD."""
    D = b * b - 4 * a * c
    P = small_primes(pmax)
    s = 1.0
    for p in P[P <= 2000]:
        p = int(p)
        w = sum(1 for n in range(p) if (a * n * n + b * n + c) % p == 0)
        s *= (1 - w / p) / (1 - 1 / p)
    big = P[P > 2000].astype(np.int64)
    d = np.mod(D, big)
    e = (big - 1) // 2
    r = np.ones_like(big)
    base = d.copy()
    while np.any(e > 0):
        m = (e & 1) == 1
        r[m] = (r[m] * base[m]) % big[m]
        base = (base * base) % big
        e >>= 1
    leg = np.where(d == 0, 0, np.where(r == 1, 1, -1))
    w = 1 + leg
    w[(a % big) == 0] = 1          # never happens for a = 4 and p > 2000
    s *= float(np.exp(np.sum(np.log1p(-w / big) - np.log1p(-1.0 / big))))
    return s
