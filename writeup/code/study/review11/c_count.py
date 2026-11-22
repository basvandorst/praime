"""Review 11, script C: own odd-only segmented sieve of [1e9, 2e9); counts of gaps 120..298 between consecutive primes
(both in the window), total number of gaps, li normalisation, band z against r48's law values."""
import math, json
import numpy as np
import mpmath

L0, L1, SEG = 10 ** 9, 2 * 10 ** 9, 10 ** 8
base = np.arange(3, math.isqrt(L1) + 2, 2)
isp = np.ones(math.isqrt(L1) + 2, bool); isp[:2] = False
for i in range(2, int(len(isp) ** .5) + 1):
    if isp[i]: isp[i * i::i] = False
base = [int(p) for p in np.nonzero(isp)[0] if p > 2]
cnt = np.zeros(400, np.int64); ngaps = 0; nprimes = 0; last = None
for lo in range(L0, L1, SEG):
    hi = lo + SEG                                    # odd numbers lo+1, lo+3, ..., hi-1 (lo even)
    odd = np.ones(SEG // 2, bool)                    # index i <-> lo + 1 + 2i
    for p in base:
        if p * p >= hi: break
        m = max(p * p, ((lo + 1 + p - 1) // p) * p)
        if m % 2 == 0: m += p
        odd[(m - lo - 1) // 2::p] = False
    P = lo + 1 + 2 * np.nonzero(odd)[0].astype(np.int64)
    nprimes += len(P)
    if last is not None: P = np.concatenate([[last], P])
    g = np.diff(P); ngaps += len(g); last = int(P[-1])
    big = g[g >= 120]
    cnt += np.bincount(big, minlength=400)[:400]
    print('segment', lo, 'primes so far', nprimes, flush=True)
li = float(mpmath.li(L1) - mpmath.li(L0)); norm = ngaps / li
print('primes', nprimes, 'gaps', ngaps, 'li diff %.3f' % li, 'norm %.10f' % norm)
r48 = json.load(open('/home/bas/praime/writeup/data/study/r48.json'))
mism = [g for g in range(120, 299, 2) if cnt[g] != r48['per_gap'][str(g)]['count']]
print('per-gap count mismatches vs r48:', mism, 'odd gaps >=120:', int(cnt[121::2].sum()), 'gaps>=300:', int(cnt[300:].sum()))
for lo, hi in ((120, 178), (180, 238), (240, 298), (120, 298)):
    c = int(sum(cnt[g] for g in range(lo, hi + 1, 2)))
    for k in ('law_o2', 'law_o3'):
        lw = sum(r48['per_gap'][str(g)][k] for g in range(lo, hi + 1, 2))
        print('band %d-%d %s count %d law %.1f ratio %.5f z_poisson %.2f' % (lo, hi, k, c, lw, c / lw, (c - lw) / math.sqrt(lw)))
    # per-gap chi-square of the band against o3
    chi = sum((cnt[g] - r48['per_gap'][str(g)]['law_o3']) ** 2 / r48['per_gap'][str(g)]['law_o3'] for g in range(lo, hi + 1, 2))
    print('   per-gap chi2 (o3) %.1f on %d gaps' % (chi, (hi - lo) // 2 + 1))
json.dump({'gaps': ngaps, 'primes': nprimes, 'norm': norm, 'cnt': {g: int(cnt[g]) for g in range(120, 400) if cnt[g]}}, open('/home/bas/scratch/study-review11/c_count.json', 'w'))
