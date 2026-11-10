# Puzzle 550: S_i = p_i p_(i+1) + p_(i+2) p_(i+3) for consecutive primes; runs of L consecutive i with S_i a semiprime.
# S_i is even, so S_i must be 2 * prime.  That needs p_i p_(i+1) = p_(i+2) p_(i+3) modulo 12 (else 3 | S_i or 4 | S_i).
# usage: p550_search.py LMIN X0 X1   (env PART=i/n splits the segments over cores).  Prints every run of >= LMIN.
import sys, os, numpy as np, gmpy2
LMIN = int(sys.argv[1]); X0 = int(float(sys.argv[2])); X1 = int(float(sys.argv[3]))
part, nparts = map(int, os.environ.get("PART", "0/1").split("/"))
SEG = 1 << 25; PAD = 20000
def small_primes(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i * i::i] = False
    return np.flatnonzero(s)
BASE = small_primes(int((X1 + SEG + PAD) ** .5) + 2)[1:]             # odd primes
def primes_between(lo, hi):                                          # lo odd; odd-only sieve
    n = (hi - lo + 1) // 2; w = np.ones(n, dtype=bool)
    for p in BASE:
        p = int(p)
        if p * p >= hi: break
        st = max(p * p, (lo + p - 1) // p * p)
        if st % 2 == 0: st += p
        w[(st - lo) // 2::p] = False
    return 2 * np.flatnonzero(w) + lo
def S(P, i): return int(P[i]) * int(P[i + 1]) + int(P[i + 2]) * int(P[i + 3])
def good(P, i):
    s = S(P, i); return s % 4 == 2 and gmpy2.is_prime(s // 2)
counts = {}; nseg = 0
for si, lo in enumerate(range(max(X0, 5) | 1, X1, SEG)):
    if si % nparts != part: continue
    P = primes_between(lo, min(lo + SEG, X1) + PAD)
    nstart = int(np.searchsorted(P, min(lo + SEG, X1)))              # starts p_i in [lo, lo + SEG)
    m = (P % 12) * (np.roll(P, -1) % 12) % 12                        # p_i p_(i+1) mod 12
    ok = m == np.roll(m, -2)
    run = ok.copy()
    for j in range(1, LMIN): run &= np.roll(ok, -j)
    idx = np.flatnonzero(run[:nstart])
    for q in (5, 7, 11, 13, 17, 19, 23, 29, 31):
        r = P % q
        for j in range(LMIN):
            if len(idx) == 0: break
            s = (r[idx + j] * r[idx + j + 1] + r[idx + j + 2] * r[idx + j + 3]) % q
            idx = idx[s != 0]
    for i in idx:
        i = int(i)
        if all(good(P, i + j) for j in range(LMIN)):
            if i > 0 and good(P, i - 1): continue                    # not the start of its run (counted from the true start)
            L = LMIN
            while good(P, i + L): L += 1
            counts[L] = counts.get(L, 0) + 1
            print("RUN", L, int(P[i]), flush=True)
    nseg += 1
    if nseg % 200 == 0: print(f"# up to {min(lo + SEG, X1):.4e}: {dict(sorted(counts.items()))}", flush=True)
print(f"# done [{X0}, {X1}) part {part}/{nparts}: runs by exact length {dict(sorted(counts.items()))}")
