#!/usr/bin/env python3
"""Exact prime sums of the Wagstaff weights, segmented odd-only numpy sieve.

  w_W(q) = e^gamma * log2(a q) / q,  a = 6 if q = 1 mod 4, a = 2 if q = 3 mod 4 (and q = 2)
  w_P(q) = e^gamma * log2(q) / q      (a = 1: the pure 'Poisson in log log' weight)

Output (npz):
  small_q, small_SW, small_SP : every prime q <= 10^6 and S(q) = sum_{r <= q} w(r)  (inclusive)
  grid_x, grid_SW, grid_SP, grid_pi : x = 1000 k (k = 1..N/1000), S(x), pi(x)
  exp_p, exp_SW, exp_SP, exp_wW, exp_wP : at each known exponent <= N
Usage: sieve_sum.py N out.npz
"""
import sys, math, time
import numpy as np

EG = math.exp(0.5772156649015329)
N = int(float(sys.argv[1]))
OUT = sys.argv[2]
EXPS = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253,
        4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243, 110503, 132049, 216091,
        756839, 859433, 1257787, 1398269, 2976221, 3021377, 6972593, 13466917, 20996011,
        24036583, 25964951, 30402457, 32582657, 37156667, 42643801, 43112609, 57885161,
        74207281, 77232917, 82589933, 136279841]
GRID = 1000
SMALL = 10**6
assert N % GRID == 0 and N >= SMALL

# base primes up to sqrt(N)
r = int(math.isqrt(N)) + 1
b = np.ones(r + 1, dtype=bool); b[:2] = False
for i in range(2, int(math.isqrt(r)) + 1):
    if b[i]:
        b[i * i::i] = False
base = np.flatnonzero(b)
base_odd = base[base > 2].astype(np.int64)

def weights(q):
    qf = q.astype(np.float64)
    a = np.where(q % 4 == 1, 6.0, 2.0)
    wW = EG * np.log2(a * qf) / qf
    wP = EG * np.log2(qf) / qf
    return wW, wP

SEGODD = 1 << 22          # odd numbers per segment -> range 2^23
totW = [EG * math.log2(2 * 2) / 2]   # q = 2 (a = 2); cancels in all differences
totP = [EG * math.log2(2) / 2]
cnt = 1
runW = math.fsum(totW); runP = math.fsum(totP)

small_q = [np.array([2], dtype=np.int64)]
small_SW = [np.array([runW])]; small_SP = [np.array([runP])]
gx = np.arange(GRID, N + 1, GRID, dtype=np.int64)
gSW = np.zeros(gx.size); gSP = np.zeros(gx.size); gpi = np.zeros(gx.size, dtype=np.int64)
exps = np.array([e for e in EXPS if e <= N], dtype=np.int64)
eSW = np.zeros(exps.size); eSP = np.zeros(exps.size); ewW = np.zeros(exps.size); ewP = np.zeros(exps.size)
# exponent 2
eSW[0] = runW; eSP[0] = runP; ewW[0] = totW[0]; ewP[0] = totP[0]

t0 = time.time()
lo = 3  # odd
# grid points below 3: x=... none (GRID=1000)
gi = 0  # next grid index to fill
while lo <= N:
    hi = min(lo + 2 * SEGODD, N + 1)      # numbers in [lo, hi)
    n = (hi - lo + 1) // 2
    seg = np.ones(n, dtype=bool)
    lim = math.isqrt(hi - 1)
    ps = base_odd[base_odd <= lim]
    if ps.size:
        st = np.maximum(ps * ps, ((lo + ps - 1) // ps) * ps)
        st = np.where(st % 2 == 0, st + ps, st)
        idx = (st - lo) // 2
        for p, i0 in zip(ps.tolist(), idx.tolist()):
            if i0 < n:
                seg[i0::p] = False
    q = lo + 2 * np.flatnonzero(seg).astype(np.int64)
    q = q[q < hi]
    wW, wP = weights(q)
    cW = np.cumsum(wW) + runW
    cP = np.cumsum(wP) + runP
    # grid points in [lo, hi): S(x) = sum over q <= x
    gj = gi
    while gj < gx.size and gx[gj] < hi:
        gj += 1
    if gj > gi:
        xs = gx[gi:gj]
        k = np.searchsorted(q, xs, side='right')
        gSW[gi:gj] = np.where(k > 0, cW[np.maximum(k - 1, 0)], runW)
        gSP[gi:gj] = np.where(k > 0, cP[np.maximum(k - 1, 0)], runP)
        gpi[gi:gj] = cnt + k
        gi = gj
    # exponents
    for j, e in enumerate(exps.tolist()):
        if lo <= e < hi:
            k = np.searchsorted(q, e)
            assert k < q.size and q[k] == e, e
            eSW[j] = cW[k]; eSP[j] = cP[k]; ewW[j] = wW[k]; ewP[j] = wP[k]
    if lo <= SMALL:
        m = q <= SMALL
        small_q.append(q[m]); small_SW.append(cW[m]); small_SP.append(cP[m])
    runW = runW + float(np.sum(wW))
    runP = runP + float(np.sum(wP))
    cnt += q.size
    lo = hi if hi % 2 == 1 else hi + 1
    if (lo // (2 * SEGODD)) % 16 == 0:
        print(f"{lo:.3e} pi={cnt} SW={runW:.9f} t={time.time()-t0:.1f}s", flush=True)

print(f"N={N} pi(N)={cnt} SW(N)={runW:.12f} SP(N)={runP:.12f} time={time.time()-t0:.1f}s")
np.savez(OUT, small_q=np.concatenate(small_q), small_SW=np.concatenate(small_SW),
         small_SP=np.concatenate(small_SP), grid_x=gx, grid_SW=gSW, grid_SP=gSP, grid_pi=gpi,
         exp_p=exps, exp_SW=eSW, exp_SP=eSP, exp_wW=ewW, exp_wP=ewP, N=N)
for x, v in [(10**7, 664579), (10**8, 5761455), (10**9, 50847534)]:
    if x <= N:
        print(f"check pi({x}) = {gpi[x // GRID - 1]} (known {v})")
