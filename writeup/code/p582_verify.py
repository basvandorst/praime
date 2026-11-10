# Puzzle 582 (follow-up question): independent check with sympy and plain integers.
# (1) Every m: m - k^2 is a square, m is composite, every prime factor p (sympy.factorint; proved below 2^64, else
#     PARI isprime) has p - k^2 a square.
# (2) Direct search over a (no products): for k = 1..20 and a <= AMAX the smallest such m, by dividing a^2 + k^2 by
#     the primes c^2 + k^2 <= a^2 + k^2.  Reproduces the page (k = 1..6) and all our values with a <= AMAX.
# (3) Two prime factors: p*q = a^2 + k^2 forces c1*c2 = k(k+1) or k(k-1) (see the write-up); all cases listed.
# Run from the repo root: (ulimit -v 1500000; .venv/bin/python3 writeup/code/p582_verify.py)   about 3 minutes
import math, numpy as np
from sympy import factorint, isprime
import cypari2
pari = cypari2.Pari()
ours = {7: 2616657418865507043448949, 8: 35674738337880689, 9: 36017834208277, 10: 3308861, 11: 2005785917,
        12: 177498661178890193089519853659369, 13: 521814892919993, 14: 34916477, 15: 102374149, 16: 152745137,
        17: 26654654355051313, 18: 131247298383974880679359973, 19: 19041556374393677, 20: 232227521}
page = {1: 10, 2: 125, 3: 5224109293, 4: 28577, 5: 334109, 6: 137677}
for k, m in sorted({**page, **ours}.items()):
    a = math.isqrt(m - k * k); assert a * a + k * k == m
    f = factorint(m); assert sum(f.values()) >= 2
    cs = []
    for p, e in f.items():
        c = math.isqrt(p - k * k); assert c * c + k * k == p and (p < 2 ** 64 and isprime(p) or pari.isprime(p)); cs += [c] * e
    print(f"k = {k}: m = {m} = {a}^2 + {k}^2 = " + " * ".join(f"({c}^2 + {k}^2)" for c in cs))

AMAX = 2 * 10 ** 6
for k in range(1, 21):
    a = np.arange(1, AMAX + 1, dtype=np.int64); m = a * a + k * k; rem = m.copy()
    good = [c * c + k * k for c in range(1, AMAX + 1) if isprime(c * c + k * k)]       # every good prime <= AMAX^2 + k^2
    for p in good:
        if p * p > AMAX * AMAX + k * k: break
        while True:
            d = rem % p == 0
            if not d.any(): break
            rem[d] //= p
    # rem is 1, or one prime factor > sqrt(m) (good only if rem - k^2 is a square), or a product of non-good primes
    r = np.sqrt((rem - k * k).clip(0).astype(np.float64)).astype(np.int64)
    cand = np.flatnonzero((rem == 1) | ((r * r + k * k == rem) & (rem < m)))
    first = None
    for i in cand:
        mm = int(m[i]); f = factorint(mm)
        if sum(f.values()) >= 2 and all(p > k * k and math.isqrt(p - k * k) ** 2 + k * k == p for p in f): first = mm; break
    known = {**page, **ours}[k]
    if known <= AMAX * AMAX: assert first == known, (k, first, known)
    else: assert first is None, (k, first)
    print(f"direct search k = {k}, a <= {AMAX}: smallest m = {first}")

for k in range(2, 21):
    sols = [(c1, n // c1) for n in (k * (k - 1), k * (k + 1)) for c1 in range(1, math.isqrt(n) + 1)
            if n % c1 == 0 and isprime(c1 * c1 + k * k) and isprime((n // c1) ** 2 + k * k)]
    assert not sols, (k, sols)
print("two prime factors: no case for k = 2..20")
print("all checks passed")
