"""Windows of integers, the primes in them, and fake primes on the same windows.

Every generator returns a boolean mask M of length W with M[i] = True when L + i is a member.
Segments of SEG numbers keep the memory low (a float array of SEG numbers is 80 MB).

Fake primes (the null ladder), all with density 1/ln n so that they are detrended like the primes:
  cramer      n is a member with probability 1/ln n (Cramer 1936).
  sieveZ      n has no prime factor <= Z, and then is a member with probability (1/ln n)/prod_{p<=Z}(1-1/p).
  class       Banks-Ford-Tao: for every prime q a class a_q mod q is struck (a_q = 0 for q <= 31, drawn by lot
              for larger q), q takes part from n_q = exp(1/prod_{q'<q}(1-1/q')) on; density ~ 1/ln n by Mertens.
              Same model as MERSENNE.md section 14 (mersenne/strategies/ten/model/ladder.py, m3 with z0 = 31).
  shuffle     the gaps of the true primes in the window, in random order (same gaps, no order).
"""
import math
import numpy as np

SEG = 10_000_000
WINDOWS = {                      # name: (L, W); holdouts are never used to choose parameters
    'train': (10_000_000, 90_000_000),
    'A': (10_000_000_000, 50_000_000),
    'B': (1_000_000_000_000, 50_000_000),
}


def small_primes(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0]


def _segments(L, W):
    lo = L
    while lo < L + W:
        hi = min(lo + SEG, L + W)
        yield lo, hi
        lo = hi


def primes_mask(L, W):
    """Exact sieve of [L, L+W)."""
    base = small_primes(int(math.isqrt(L + W)) + 1)
    M = np.ones(W, dtype=bool)
    for lo, hi in _segments(L, W):
        seg = M[lo - L:hi - L]
        for p in base:
            p = int(p)
            if p * p >= hi:
                break
            st = max(p * p, -(-lo // p) * p)
            seg[st - lo::p] = False
    if L < 2:
        M[:2 - L] = False
    return M


def _lnvals(lo, hi):
    return np.log(np.arange(lo, hi, dtype=np.float64))


def cramer_mask(L, W, rng):
    M = np.empty(W, dtype=bool)
    for lo, hi in _segments(L, W):
        M[lo - L:hi - L] = rng.random(hi - lo) < 1.0 / _lnvals(lo, hi)
    return M


def sieve_mask(L, W, Z, rng):
    P = small_primes(Z)
    boost = 1.0 / np.prod([1.0 - 1.0 / p for p in P])
    M = np.ones(W, dtype=bool)
    for p in P:
        p = int(p)
        M[(-L) % p::p] = False
    for lo, hi in _segments(L, W):
        seg = M[lo - L:hi - L]
        seg &= rng.random(hi - lo) < boost / _lnvals(lo, hi)
    return M


def class_mask(L, W, rng, z0=31):
    hi = L + W
    # largest q needed: n_q <= hi  <=>  prod_{q'<q}(1-1/q') >= 1/ln hi ; Mertens gives q ~ hi^0.5615
    zguess = int(math.exp(0.5615 * math.log(hi)) * 1.3) + 100
    Q = small_primes(zguess)
    M = np.ones(W, dtype=bool)
    P = 1.0
    for q in Q:
        q = int(q)
        nq = math.exp(1.0 / P) if P < 1 else 2.0
        if nq > hi:
            break
        P *= 1.0 - 1.0 / q
        start = max(L, int(math.ceil(nq)))
        a = 0 if q <= z0 else int(rng.integers(q))
        st = start + ((a - start) % q)
        if q <= z0 and st == q:
            st += q
        M[st - L::q] = False
    else:
        raise RuntimeError('zguess too small')
    if L < 2:
        M[:2 - L] = False
    return M


def shuffle_mask(Mp, L, rng):
    P = np.nonzero(Mp)[0]
    g = np.diff(P)
    rng.shuffle(g)
    Q = np.concatenate(([P[0]], P[0] + np.cumsum(g)))
    M = np.zeros_like(Mp)
    M[Q] = True
    return M


GENS = ['primes', 'cramer', 'sieve13', 'sieve47', 'sieve100', 'class', 'shuffle']


def make_mask(gen, L, W, rng, Mp=None):
    if gen == 'primes':
        return primes_mask(L, W)
    if gen == 'cramer':
        return cramer_mask(L, W, rng)
    if gen.startswith('sieve'):
        return sieve_mask(L, W, int(gen[5:]), rng)
    if gen == 'class':
        return class_mask(L, W, rng)
    if gen == 'shuffle':
        return shuffle_mask(Mp if Mp is not None else primes_mask(L, W), L, rng)
    raise ValueError(gen)


def li(x):
    """Logarithmic integral li(x) for x > 1 (Ramanujan's series)."""
    x = float(x)
    lx = math.log(x)
    s, term, inner = 0.0, 1.0, 0.0
    for n in range(1, 200):
        term *= lx / n
        if (n - 1) % 2 == 0:
            inner = sum(1.0 / (2 * k + 1) for k in range(0, (n - 1) // 2 + 1))
        s += ((-1) ** (n - 1)) * term / (2 ** (n - 1)) * inner
    return 0.5772156649015329 + math.log(lx) + math.sqrt(x) * s
