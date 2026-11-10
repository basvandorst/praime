# Puzzles 655 / 498, calibration of the density model by plain enumeration.
# All consecutive primes p < q < r < s <= X (segmented numpy sieve, 400 MB): N = pq + rs (type 1), ps + qr (2),
# pr + qs (3).  Lists the exact squares and counts the near-squares N = m^2 + 5040 j, |j| <= 500 (same residues
# mod 16, 9, 5, 7 as a square), for p >= 2e6.  The model: N is at a given distance from a square with probability
# 1/(2 sqrt N) = 1/(2 sqrt 2 p); the ratio counted/expected measures its error.
# Usage: p655_calib.py X      (X = 3e9: 20 s)
import sys, math, numpy as np
X = int(float(sys.argv[1])); SEG = 1 << 25; Jm = 500; P0 = 2_000_000
base = np.ones(math.isqrt(X) + 2, dtype=bool); base[:2] = False
for i in range(2, math.isqrt(len(base)) + 1):
    if base[i]: base[i * i::i] = False
base = np.flatnonzero(base)
carry = np.array([], dtype=np.int64)
exact = {1: [], 2: [], 3: []}; near = {1: 0, 2: 0, 3: 0}; naive = 0.0; hist = {1: {}, 2: {}, 3: {}}
for lo in range(0, X, SEG):
    n = min(SEG, X - lo); w = np.ones(n, dtype=bool)
    for p in base:
        p = int(p)
        if p * p >= lo + n: break
        w[max(p * p, (lo + p - 1) // p * p) - lo::p] = False
    if lo == 0: w[:2] = False
    P = np.concatenate((carry, np.flatnonzero(w) + lo)); carry = P[-3:]
    if len(P) < 4: continue
    p, q, r, s = (P[i:len(P) - 3 + i] for i in range(4))
    pu, qu, ru, su = (x.astype(np.uint64) for x in (p, q, r, s)); pf, qf, rf, sf = (x.astype(np.float64) for x in (p, q, r, s))
    big = p >= P0
    naive += float((1.0 / (2.0 * np.sqrt(2.0) * pf[big])).sum())
    for typ, (Nu, Nf) in {1: (pu * qu + ru * su, pf * qf + rf * sf), 2: (pu * su + qu * ru, pf * sf + qf * rf), 3: (pu * ru + qu * su, pf * rf + qf * sf)}.items():
        m = np.rint(np.sqrt(Nf)).astype(np.uint64)
        delta = (Nu - m * m).astype(np.int64)
        for i in np.flatnonzero(delta == 0):
            a, b, c, d = (int(x[i]) for x in (p, q, r, s)); N = {1: a * b + c * d, 2: a * d + b * c, 3: a * c + b * d}[typ]
            if math.isqrt(N) ** 2 == N: exact[typ].append((a, b, c, d))
        sel = big & (delta % 5040 == 0) & (np.abs(delta) <= 5040 * Jm)
        near[typ] += int(sel.sum())
print("X =", X, "exact:", exact)
import math as _m
print(f"model: expected exact squares per type for 2e6 <= p <= X: {naive:.3f}; 0.354 ln(ln X / ln 2e6) = {0.3536 * _m.log(_m.log(X) / _m.log(2e6)):.3f}")
print("near-squares (p >= %d): %s; naive expectation %.1f each; ratio %s" % (P0, near, naive * (2 * Jm + 1), {t: round(near[t] / (naive * (2 * Jm + 1)), 3) for t in near}))
