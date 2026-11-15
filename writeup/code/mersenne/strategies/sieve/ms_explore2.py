"""Exploration: is P(2 is a 2k-th power residue mod q) = 1/k for primes q = 1+2k mod 4k, q = +-1 mod 8?
Population A: all such primes q < 4e7.  Population B: those with (q-1)/(2k) prime (q = 2kp+1).
Vectorised modpow in numpy (q < 2^26, products < 2^52)."""
import numpy as np, math
X = 40_000_000
s = np.ones(X + 1, dtype=bool); s[:2] = False
for i in range(2, int(X ** 0.5) + 1):
    if s[i]: s[i * i::i] = False
Q = np.nonzero(s)[0].astype(np.uint64)
def powmod(base, e, m):
    r = np.ones(len(m), dtype=np.uint64); b = np.full(len(m), base, dtype=np.uint64) % m
    e = e.copy()
    while np.any(e):
        odd = (e & np.uint64(1)).astype(bool)
        r = np.where(odd, (r * b) % m, r)
        b = (b * b) % m
        e >>= np.uint64(1)
    return r
tot = {}
for k in [3, 4, 5, 7, 8, 9, 11, 12, 13, 15, 16, 17, 20, 21, 24, 28, 32, 36, 40, 48, 60, 64]:
    sel = (Q % np.uint64(4 * k) == np.uint64(1 + 2 * k)) & ((Q % np.uint64(8) == 1) | (Q % np.uint64(8) == 7))
    q = Q[sel]
    hit = powmod(2, (q - np.uint64(1)) // np.uint64(2 * k), q) == 1
    pB = (q - np.uint64(1)) // np.uint64(2 * k)
    isB = s[pB.astype(np.int64)]
    nA, hA = len(q), int(hit.sum()); nB, hB = int(isB.sum()), int(hit[isB].sum())
    print(f"k={k:3d}: all q: {hA}/{nA} = {hA/nA:.4f} (1/k={1/k:.4f}, z={(hA-nA/k)/math.sqrt(nA*(1/k)*(1-1/k)):+.2f})"
          f"   q=2kp+1: {hB}/{nB} = {hB/nB:.4f} (z={(hB-nB/k)/math.sqrt(nB*(1/k)*(1-1/k)):+.2f})")
