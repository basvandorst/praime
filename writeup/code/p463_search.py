# Puzzle 463: T_k = k p_k + (k+1) p_(k+1) (p_k the k-th prime); runs of consecutive k with T_k prime.
# Segmented sieve with the index k carried along (pi(X0) from sympy.primepi); T_k is tested modulo 3..31 on residues first.
# usage: p463_search.py LMIN X0 X1      prints every run of >= LMIN as "RUN length k p_k"
import sys, numpy as np, gmpy2
from sympy import primepi
LMIN = int(sys.argv[1]); X0 = int(float(sys.argv[2])); X1 = int(float(sys.argv[3]))
SEG = 1 << 25; PAD = 20000
def small_primes(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i * i::i] = False
    return np.flatnonzero(s)
BASE = small_primes(int((X1 + SEG + PAD) ** .5) + 2)[1:]
def primes_between(lo, hi):                                          # lo odd, lo >= 3
    n = (hi - lo + 1) // 2; w = np.ones(n, dtype=bool)
    for p in BASE:
        p = int(p)
        if p * p >= hi: break
        st = max(p * p, (lo + p - 1) // p * p)
        if st % 2 == 0: st += p
        w[(st - lo) // 2::p] = False
    return 2 * np.flatnonzero(w) + lo
def T(k, P, i): return k * int(P[i]) + (k + 1) * int(P[i + 1])
lo0 = max(X0, 3) | 1
index = int(primepi(lo0 - 1))                                       # number of primes below the first segment
counts = {}; nseg = 0
for lo in range(lo0, X1, SEG):
    P = primes_between(lo, min(lo + SEG, X1) + PAD)
    nstart = int(np.searchsorted(P, min(lo + SEG, X1)))
    K = index + 1 + np.arange(len(P), dtype=np.int64)                # P[i] is the K[i]-th prime
    ok = np.ones(len(P), dtype=bool)
    for q in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
        r = P % q; kk = K % q
        ok &= (kk * r + (kk + 1) * np.roll(r, -1)) % q != 0
    run = ok.copy()
    for j in range(1, LMIN): run &= np.roll(ok, -j)
    for i in np.flatnonzero(run[:nstart]):
        i = int(i); k = int(K[i])
        if all(gmpy2.is_prime(T(k + j, P, i + j)) for j in range(LMIN)):
            if i > 0 and gmpy2.is_prime(T(k - 1, P, i - 1)): continue   # counted from its true start
            L = LMIN
            while gmpy2.is_prime(T(k + L, P, i + L)): L += 1
            counts[L] = counts.get(L, 0) + 1
            print("RUN", L, k, int(P[i]), flush=True)
    index += nstart; nseg += 1
    if nseg % 200 == 0: print(f"# up to {min(lo + SEG, X1):.4e}: {dict(sorted(counts.items()))}", flush=True)
print(f"# done [{X0}, {X1}): runs by exact length {dict(sorted(counts.items()))}")
