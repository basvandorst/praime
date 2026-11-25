"""Ideals of Z[i] with norm <= N, with angle, content, Moebius, divisor count, primality.

Each nonzero ideal is represented once by its generator a + bi with a >= 1, b >= 0.
Written independently for the referee report (no code taken from rh_check.py).
"""
import math
import numpy as np


def _gauss_divides(a, b, p, q):
    """True if (p + qi) divides (a + bi) in Z[i]."""
    n = p * p + q * q
    # (a+bi)(p-qi) = (ap+bq) + (bp-aq)i
    re = a * p + b * q
    im = b * p - a * q
    return re % n == 0 and im % n == 0


def _gauss_div(a, b, p, q):
    n = p * p + q * q
    re = a * p + b * q
    im = b * p - a * q
    return re // n, im // n


def _rational_factor(n):
    f = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


_split_cache = {}


def _split_prime(p):
    """A Gaussian prime u + vi of norm p (p = 2 or p = 1 mod 4)."""
    if p in _split_cache:
        return _split_cache[p]
    for u in range(1, int(math.isqrt(p)) + 1):
        v2 = p - u * u
        v = math.isqrt(v2)
        if v * v == v2:
            _split_cache[p] = (u, v)
            return (u, v)
    raise ValueError(p)


def ideal_factor(a, b):
    """Exponents of prime ideals in (a + bi): list of (key, exponent)."""
    n = a * a + b * b
    out = []
    for p, e in _rational_factor(n).items():
        if p == 2:
            out.append((('2',), e))  # (1+i)^e, N(1+i) = 2
        elif p % 4 == 3:
            out.append((('q', p), e // 2))
        else:
            u, v = _split_prime(p)
            # count powers of pi = u + vi in (a+bi)
            x, y, k = a, b, 0
            while _gauss_divides(x, y, u, v):
                x, y = _gauss_div(x, y, u, v)
                k += 1
            if k:
                out.append((('p', p, 0), k))
            if e - k:
                out.append((('p', p, 1), e - k))
    return out


def enumerate_ideals(N, lo=0):
    """Ideals with lo < N(a) <= N. Returns dict of numpy arrays."""
    A, B = [], []
    r = math.isqrt(N)
    for a in range(1, r + 1):
        for b in range(0, r + 1):
            n = a * a + b * b
            if lo < n <= N:
                A.append(a)
                B.append(b)
    A = np.array(A, dtype=np.int64)
    B = np.array(B, dtype=np.int64)
    norm = A * A + B * B
    theta = np.arctan2(B.astype(float), A.astype(float))  # in [0, pi/2)
    cont = np.gcd(A, B)
    prim = (A // cont) * 100003 + (B // cont)  # id of primitive part
    mu = np.zeros(len(A), dtype=np.int64)
    dv = np.zeros(len(A), dtype=np.int64)
    isprime = np.zeros(len(A), dtype=bool)
    for i in range(len(A)):
        fac = ideal_factor(int(A[i]), int(B[i]))
        d = 1
        sq = False
        for _, e in fac:
            d *= e + 1
            if e >= 2:
                sq = True
        dv[i] = d
        mu[i] = 0 if sq else (-1) ** len(fac)
        isprime[i] = (len(fac) == 1 and fac[0][1] == 1)
    return dict(a=A, b=B, norm=norm, theta=theta, cont=cont, prim=prim, mu=mu, d=dv,
                isprime=isprime)
