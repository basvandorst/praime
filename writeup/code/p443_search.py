# Puzzles 443 and 393 (Q1): is a sum of k-th powers of L >= 2 consecutive primes ever a k-th power (k = 3)?
# Every window sum N <= X is tested.  The sums are built by the recurrence S_L(i) = S_(L-1)(i) + p_(i+L-1)^k, once
# in float64 (all terms positive, so no cancellation: relative error below L * 2^-53) and once modulo 2^64.
# m = round(N^(1/k)) is exact as long as the float error is below k/2 * N^(-1/k) relative; then N = m^k iff
# m^k = N (mod 2^64), up to a chance of 2^-64, and every hit is recomputed with Python integers.
# The script also adds up the heuristic expectation: sum over the windows of (1/k) N^(1/k - 1), per decade of N.
# Usage: p443_search.py k X      (ulimit -v 1500000; memory about 150 MB for k = 3, X = 1e24)
import sys, math
import numpy as np
from sympy import integer_nthroot
k = int(sys.argv[1]); X = float(sys.argv[2])
PMAX = int((X / 2) ** (1.0 / k)) + 100                       # a window of 2 primes has both below (X/2)^(1/k) + gap
sv = np.ones(PMAX + 1, dtype=bool); sv[:2] = False
for i in range(2, math.isqrt(PMAX) + 1):
    if sv[i]: sv[i * i::i] = False
P = np.flatnonzero(sv); del sv
Pf = P.astype(np.float64) ** k
Pw = P.astype(np.uint64) ** np.uint64(k)                     # wraps modulo 2^64
Sf = Pf.copy(); Sw = Pw.copy()                               # windows of length 1
hits = []; tested = 0; maxL = 1
dec = {}                                                     # decade of N -> [windows, expected]
L = 1
while True:
    L += 1
    n = len(Sf) - 1
    if n <= 0: break
    Sf = Sf[:n] + Pf[L - 1:L - 1 + n]
    Sw = Sw[:n] + Pw[L - 1:L - 1 + n]
    n = int(np.searchsorted(Sf, X, side='right'))            # Sf increases with the start
    if n == 0: break
    Sf = Sf[:n]; Sw = Sw[:n]; tested += n; maxL = L
    m = np.rint(Sf ** (1.0 / k)) if k != 3 else np.rint(np.cbrt(Sf))
    mw = m.astype(np.uint64)
    for d in (-1, 0, 1):                                     # +-1 guards the rounding of the float root
        md = mw + np.uint64(d) if d >= 0 else mw - np.uint64(1)
        for i in np.flatnonzero(md ** np.uint64(k) == Sw):
            N = sum(int(p) ** k for p in P[i:i + L])
            r, exact = integer_nthroot(N, k)
            if exact: hits.append((r, L, int(P[i]), int(P[i + L - 1])))
    e = Sf ** (1.0 / k - 1.0) / k
    lg = np.floor(np.log10(Sf)).astype(np.int64)
    for dd in np.unique(lg):
        sel = lg == dd
        t = dec.setdefault(int(dd), [0, 0.0, 0.0])
        t[0] += int(sel.sum()); t[1] += float(e[sel].sum())
        if L >= 3: t[2] += float(e[sel].sum())
assert maxL * 2.0 ** -53 < 0.1 * k / 2 * X ** (-1.0 / k), "float sums too coarse for exact roots"
print(f"k = {k}, X = {X:.0e}: {len(P)} primes below {PMAX}, {tested} window sums of 2..{maxL} consecutive primes tested")
tot = 0.0; tot3 = 0.0
for d in sorted(dec):
    tot += dec[d][1]; tot3 += dec[d][2]
print(f"heuristic expectation of solutions with N <= X: {tot:.3f} (L >= 3 only: {tot3:.3f})")
ds = sorted(dec)[-6:]
print("last decades, 10^d <= N < 10^(d+1):  " + "  ".join(f"d={d}: {dec[d][0]} windows, {dec[d][1]:.2e}" for d in ds))
if k >= 3:
    r = (dec[ds[-1]][1] / dec[ds[-4]][1]) ** (1 / 3)         # mean ratio of the last three steps
    print(f"ratio per decade {r:.3f}; geometric tail beyond X: {dec[ds[-1]][1] * r / (1 - r):.1e}")
print(f"solutions: {len(set(hits))}")
for r, L, a, b in sorted(set(hits)): print(f"  {a}^{k} + ... + {b}^{k} ({L} primes) = {r}^{k}")
