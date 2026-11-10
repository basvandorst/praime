# Puzzle 582 (follow-up question): composite m = a^2 + k^2 whose prime divisors are all of the form c^2 + k^2.
# Lists every product m <= M of at least NMIN such primes and tests whether m - k^2 is a square
# (root in floating point, exact test modulo 2^64, then exact integers).  Valid for M < 1e27.
# usage: p582_products.py k M NMIN      e.g.  p582_products.py 8 100000000000000000 3
import sys, math, bisect, numpy as np, gmpy2
k = int(sys.argv[1]); M = int(sys.argv[2]); NMIN = int(sys.argv[3]); K2 = k * k
pmin = next(c * c + K2 for c in range(1, 10 ** 5) if gmpy2.is_prime(c * c + K2))
cmax = math.isqrt(M // pmin ** (NMIN - 1)) + 1
gl = [c * c + K2 for c in range(1, cmax + 1) if gmpy2.is_prime(c * c + K2)]
gu = np.array(gl, dtype=np.uint64); gf = gu.astype(np.float64)
MASK = (1 << 64) - 1
hits = []; leaves = 0
def rec(prod, idx, nf):
    global leaves
    if nf + 1 >= NMIN:
        hi = bisect.bisect_right(gl, M // prod)
        if hi > idx:
            leaves += hi - idx
            r = np.rint(np.sqrt(float(prod) * gf[idx:hi] - K2)).astype(np.uint64)
            ok = r * r + np.uint64(K2) == np.uint64(prod & MASK) * gu[idx:hi]
            for j in np.flatnonzero(ok):
                m = prod * gl[idx + int(j)]; a = math.isqrt(m - K2)
                if a * a + K2 == m: hits.append(m)
    for j in range(idx, len(gl)):
        p = gl[j]
        if prod * p ** max(2, NMIN - nf) > M: break
        rec(prod * p, j, nf + 1)
with np.errstate(over="ignore"): rec(1, 0, 0)
print("k", k, "M", M, "products with at least", NMIN, "prime factors:", leaves, "hits", sorted(hits))
