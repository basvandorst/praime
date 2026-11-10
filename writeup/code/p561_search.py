# Puzzle 561: runs of L consecutive primes p_0 < ... < p_(L-1) and an even k with p_i^2 + (k + 2i)^2 prime for all i.
# usage: p561_search.py L PMIN PMAX   (env PART=i/n splits the prime segments over cores)
import sys, os, math, numpy as np, gmpy2
L = int(sys.argv[1]); PMIN = int(float(sys.argv[2])); PMAX = int(float(sys.argv[3]))
part, nparts = map(int, os.environ.get("PART", "0/1").split("/"))
WHEEL = [5, 13, 17, 29, 37, 41]
Q = math.prod(WHEEL); M = 2 * Q                                   # k runs over one period 0 < k <= M, k even
def small_primes(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i * i::i] = False
    return np.flatnonzero(s)
SP = small_primes(40000)
SIEVE = [int(q) for q in SP if q % 4 == 1 and 41 < q < 1000]
ROOT = {q: next(x for x in range(2, q) if (x * x + 1) % q == 0) for q in WHEEL + SIEVE}
E = {q: (M // q) * pow(M // q, -1, q) % M for q in WHEEL}          # = 1 mod q, 0 mod the other wheel primes and mod 2
IDX = np.arange(L, dtype=np.int64)
# constant part of the Bateman-Horn product
CONST = 2.0 ** L
for q in small_primes(10 ** 7)[1:]:
    q = int(q)
    if q % 4 == 3: CONST *= (1 - 1 / q) ** -L
    elif q in WHEEL: CONST *= (1 - 1 / q) ** -L                    # the wheel has removed the bad residues
    elif q > 1000: CONST *= (1 - 2 * L / q) * (1 - 1 / q) ** -L
def excluded(q, pr):
    s = ROOT[q]; a = (s * pr - 2 * IDX) % q; b = (-s * pr - 2 * IDX) % q
    t = np.zeros(q, dtype=bool); t[a] = True; t[b] = True
    return t
def segment_primes(lo, hi):
    w = np.ones(hi - lo, dtype=bool)
    for p in SP:
        p = int(p)
        if p * p >= hi: break
        w[max(p * p, (lo + p - 1) // p * p) - lo::p] = False
    return np.flatnonzero(w) + lo
SEG = 2 * 10 ** 7
expect = 0.0; ncand = 0; nruns = 0; hits = []
carry = np.zeros(0, dtype=np.int64)
for si, lo in enumerate(range(PMIN, PMAX, SEG)):
    if si % nparts != part: continue
    P = segment_primes(lo, min(lo + SEG, PMAX) + 2000)             # overlap so that runs starting in the segment are complete
    n = int(np.searchsorted(P, min(lo + SEG, PMAX)))
    if len(P) < n + L: n = len(P) - L
    r = P % 5; mask = np.zeros(n, dtype=np.int64)
    for i in range(L):
        ri = r[i:i + n]
        mask |= (1 << ((2 * ri - 2 * i) % 5)) | (1 << ((-2 * ri - 2 * i) % 5))
    for j in np.flatnonzero(mask != 31):
        pr = P[j:j + L]
        tabs = [excluded(q, pr) for q in WHEEL]
        if any(t.all() for t in tabs): continue
        ks = np.zeros(1, dtype=np.int64)
        for q, t in zip(WHEEL, tabs):
            ks = ((ks[:, None] + (np.flatnonzero(~t) * E[q])[None, :]) % M).ravel()
        ks[ks == 0] = M
        nruns += 1; ncand += len(ks)
        c = CONST
        lnsum = (np.log(float(pr[0]) ** 2 + ks.astype(np.float64) ** 2) ** -L).sum()   # sum of 1/(ln N)^L over the wheel survivors
        for q in SIEVE:
            t = excluded(q, pr); c *= (1 - t.sum() / q) * (1 - 1 / q) ** -L
            if len(ks): ks = ks[~t[ks % q]]
        expect += c * lnsum
        for k in ks:
            k = int(k)
            if all(gmpy2.is_prime(int(p) ** 2 + (k + 2 * i) ** 2) for i, p in enumerate(pr)):
                hits.append((int(pr[0]), k)); print("HIT", L, int(pr[0]), k, flush=True)
    print(f"p < {min(lo + SEG, PMAX):.3e}: runs {nruns}, candidates {ncand:.3e}, expected {expect:.4f}, found {len(hits)}", flush=True)
