"""The perfect row: Mersenne exponents one after another from 2, each proved. For every prime p in order:
trial factors f = 2kp + 1 (k <= K(p), f = +-1 mod 8, no prime factor <= 3000; 2^p mod f by the exact float-assisted
uint64 powering of writeup/code/mersenne/touch/tfdeep.py), then Lucas-Lehmer in gmpy2 (s <- s^2 - 2 mod 2^p - 1,
reduced by shift and add). A factor proves 2^p - 1 composite; s_(p-2) = 0 proves it prime.
Prints one line per exponent: number, p, wall time, primes tested so far, LL tests so far.  usage: row.py [stop_p]"""
import sys, time, math, numpy as np, gmpy2
from sympy import nextprime

def pow2_mod(p, q):
    qi = q.view(np.int64); qf = q.astype(np.float64)
    r = np.ones(len(q), dtype=np.uint64)
    for b in bin(p)[2:]:
        h = np.floor(r.astype(np.float64) * r.astype(np.float64) / qf).astype(np.uint64)
        t = (r * r - h * q).view(np.int64)
        t = np.where(t < 0, t + qi, t); t = np.where(t >= qi, t - qi, t)
        if b == '1':
            t = t << 1; t = np.where(t >= qi, t - qi, t)
        r = t.view(np.uint64)
    return r

SP = [l for l in range(3, 3001) if all(l % d for d in range(2, int(l ** 0.5) + 1))]

def factor(p):
    K = max(1000, min(int(1e5 * (p / 7000) ** 2), 10 ** 6))
    k = np.arange(1, K + 1, dtype=np.uint64)
    q = np.uint64(2 * p) * k + np.uint64(1)
    m = (q % np.uint64(8) == 1) | (q % np.uint64(8) == 7)
    for l in SP:
        m &= (q % np.uint64(l) != 0) | (q == np.uint64(l))
    q = q[m]
    if p < 64:
        q = q[q < np.uint64((1 << p) - 1)]
    if len(q) == 0: return 0
    hit = np.flatnonzero(pow2_mod(p, q) == 1)
    return int(q[hit[0]]) if len(hit) else 0

def ll(p):
    if p == 2: return True
    M = (gmpy2.mpz(1) << p) - 1; s = gmpy2.mpz(4)
    for _ in range(p - 2):
        s = s * s - 2
        s = (s & M) + (s >> p)
        if s >= M: s -= M
    return s == 0

stop = int(sys.argv[1]) if len(sys.argv) > 1 else 10 ** 9
t0 = time.time(); p = 2; n = 0; tested = 0; lls = 0
while p <= stop:
    tested += 1
    if p == 2 or factor(p) == 0:
        lls += 1
        if ll(p):
            n += 1
            print(f'{n:2d}. {p:>8}  {time.time() - t0:9.1f} s  primes tested {tested:,}  Lucas-Lehmer tests {lls:,}', flush=True)
    p = int(nextprime(p))
print(f'stopped at p = {p}: {n} exponents, {tested:,} primes, {lls:,} Lucas-Lehmer tests, {time.time() - t0:.0f} s', flush=True)
