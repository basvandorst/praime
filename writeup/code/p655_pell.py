# Puzzles 655 and 498: consecutive primes p < q < r < s with pq + rs (type 1), ps + qr (type 2) or pr + qs (type 3)
# a perfect square.  With q = p + 2a, r = p + 2b, s = p + 2d (half-gaps a < b < d), t = a + b + d and
# C = bd, ab, ad for the three types:   N = 2p^2 + 2tp + 4C = (2z)^2   <=>   U^2 - 8 z^2 = k,  U = 2p + t, k = t^2 - 8C,
# and t must be odd.  All solutions of U^2 - 8z^2 = k are the orbits of finitely many fundamental solutions under
# (U, z) -> (3U + 8z, U + 3z).  So: list the fundamental solutions for every k, walk the orbits (5.83 per step),
# and for every (orbit point, gap triple with that k) test p = (U - t)/2.  The orbit is carried modulo
# Q = 3*5*...*23 and Q2 = 29*...*43 in numpy to sieve p, q, r, s; survivors are rebuilt as exact integers and tested with gmpy2.
# Usage: p655_pell.py NDIG [DMAX] [PMIN] [--raw]     (--raw: no sieve, no primality: all integer solutions, for tests)
# PART=i/n in the environment runs every n-th chunk of d only (to use several cores).
# Memory: chunks of at most 3e6 (orbit, triple) pairs, below 600 MB.  Run under ulimit -v 1500000.
import sys, math, time
import numpy as np
import gmpy2
args = [a for a in sys.argv[1:] if not a.startswith('--')]; RAW = '--raw' in sys.argv
NDIG = int(args[0]); DMAX = int(args[1]) if len(args) > 1 else 7 * NDIG + 30
PMIN = int(float(args[2])) if len(args) > 2 else 3
QP = [3, 5, 7, 11, 13, 17, 19, 23]; Q = math.prod(QP); INV2 = (Q + 1) // 2
QP2 = [29, 31, 37, 41, 43]; Q2 = math.prod(QP2); INV22 = (Q2 + 1) // 2      # second residue track
def dcap(digits): return 7 * digits + 30            # largest half-gap sum d used at this size: s - p <= 6 ln p + 60

# ---- fundamental solutions of U^2 - 8 z^2 = k for all odd |k| <= KMAX (Nagell's bounds: 8 z^2 <= k for k > 0,
#      4 z^2 <= |k| for k < 0), both signs of z
KMAX = 9 * DMAX * DMAX
fk, fU, fz = [], [], []
for z in range(0, math.isqrt(KMAX // 4) + 1):
    U = np.arange(1, math.isqrt(KMAX + 8 * z * z) + 2, 2, dtype=np.int64)
    k = U * U - 8 * z * z
    ok = (np.abs(k) <= KMAX) & np.where(k > 0, 8 * z * z <= k, 4 * z * z <= -k)
    for sgn in ((1,) if z == 0 else (1, -1)):
        fk.append(k[ok]); fU.append(U[ok]); fz.append(np.full(int(ok.sum()), sgn * z, dtype=np.int64))
fk, fU, fz = (np.concatenate(x) for x in (fk, fU, fz))
# A fundamental solution with z < 0 can step to U < 0; (U, z) and (-U, -z) are the same solution.  Walk such
# entries until U > 0 and z >= 0 (from there the orbit stays positive); the points passed on the way are kept as
# single points (fone) so that no value of U is lost.
sU, sz, sk, s1 = [fU[fz >= 0]], [fz[fz >= 0]], [fk[fz >= 0]], [np.zeros(int((fz >= 0).sum()), dtype=bool)]
cU, cz, ck = fU[fz < 0], fz[fz < 0], fk[fz < 0]
while len(cU):
    sU.append(cU); sz.append(cz); sk.append(ck); s1.append(np.ones(len(cU), dtype=bool))      # the point itself
    cU, cz = 3 * cU + 8 * cz, cU + 3 * cz
    neg = cU < 0; cU = np.where(neg, -cU, cU); cz = np.where(neg, -cz, cz)
    done = cz >= 0
    sU.append(cU[done]); sz.append(cz[done]); sk.append(ck[done]); s1.append(np.zeros(int(done.sum()), dtype=bool))
    cU, cz, ck = cU[~done], cz[~done], ck[~done]
fU, fz, fk, fone = (np.concatenate(x) for x in (sU, sz, sk, s1))
o = np.argsort(fk, kind='stable'); fk, fU, fz, fone = fk[o], fU[o], fz[o], fone[o]
NSTEP = int((NDIG + 1) / math.log10(3 + math.sqrt(8))) + 2
E = [(1, 0)]                                         # (3 + sqrt 8)^n = A + B sqrt 8
for n in range(NSTEP): A, B = E[-1]; E.append((3 * A + 8 * B, A + 3 * B))
LOG = math.log10(3 + math.sqrt(8))
print(f"NDIG = {NDIG}, half-gap d <= {DMAX}, |k| <= {KMAX}: {len(fk)} fundamental solutions, {NSTEP} orbit steps", flush=True)

def triples(dlo, dhi):
    """records (typ, a, b, d, t, k) for dlo <= d < dhi: t odd, {0, a, b, d} not all classes mod 3"""
    out = []
    for d in range(dlo, dhi):
        a, b = np.meshgrid(np.arange(1, d, dtype=np.int64), np.arange(1, d, dtype=np.int64), indexing='ij')
        sel = (a < b) & ((a + b + d) % 2 == 1)
        a = a[sel]; b = b[sel]
        res = [(0 == r) | ((2 * a) % 3 == r) | ((2 * b) % 3 == r) | ((2 * d) % 3 == r) for r in range(3)]
        adm = ~(res[0] & res[1] & res[2])
        if RAW: adm[:] = True
        a = a[adm]; b = b[adm]; t = a + b + d; dd = np.full(len(a), d, dtype=np.int64)
        for typ, C in ((1, b * d), (2, a * b), (3, a * d)):
            out.append((np.full(len(a), typ, dtype=np.int64), a, b, dd, t, t * t - 8 * C))
    return [np.concatenate([o[i] for o in out]) for i in range(6)] if out else None

solutions = []; stats = {}; npairs = 0; nprp = 0
def run_chunk(rec):
    global npairs, nprp
    typ, a, b, d, t, k = rec
    lo = np.searchsorted(fk, k, side='left'); hi = np.searchsorted(fk, k, side='right'); cnt = hi - lo
    keep = cnt > 0
    if not keep.any(): return
    ridx = np.repeat(np.flatnonzero(keep), cnt[keep])                      # record of each pair
    start = np.repeat(lo[keep], cnt[keep]); first = np.cumsum(cnt[keep]) - cnt[keep]
    fidx = start + (np.arange(len(ridx)) - np.repeat(first, cnt[keep]))    # fundamental solution of each pair
    used, inv = np.unique(fidx, return_inverse=True)
    Um = fU[used] % Q; zm = fz[used] % Q; Um2 = fU[used] % Q2; zm2 = fz[used] % Q2; ulog = np.log10(fU[used].astype(np.float64) + 3.0 * np.abs(fz[used]))
    pt, pa, pb, pd = t[ridx], a[ridx], b[ridx], d[ridx]
    for n in range(NSTEP + 1):
        if n:
            Um, zm = (3 * Um + 8 * zm) % Q, (Um + 3 * zm) % Q
            Um2, zm2 = (3 * Um2 + 8 * zm2) % Q2, (Um2 + 3 * zm2) % Q2
        dig = ulog + n * LOG - math.log10(2.0)                             # about log10 p
        act = (dig[inv] <= NDIG + 0.5) & (dig[inv] >= math.log10(PMIN) - 1.0) & (pd <= dcap(np.maximum(dig[inv], 1.0)))
        if n: act &= ~fone[used][inv]
        idx = np.flatnonzero(act)
        if len(idx) == 0: continue
        npairs += len(idx)
        if not RAW:
            pm = ((Um[inv[idx]] - pt[idx]) % Q) * INV2 % Q
            for q in QP:
                r = pm % q
                ok = (r != 0) & ((r + 2 * pa[idx]) % q != 0) & ((r + 2 * pb[idx]) % q != 0) & ((r + 2 * pd[idx]) % q != 0)
                idx = idx[ok]; pm = pm[ok]
                if len(idx) == 0: break
            pm = ((Um2[inv[idx]] - pt[idx]) % Q2) * INV22 % Q2
            for q in QP2:
                r = pm % q
                ok = (r != 0) & ((r + 2 * pa[idx]) % q != 0) & ((r + 2 * pb[idx]) % q != 0) & ((r + 2 * pd[idx]) % q != 0)
                idx = idx[ok]; pm = pm[ok]
        A, B = E[n]
        for i in idx:
            f = int(used[inv[i]]); U = A * int(fU[f]) + 8 * B * int(fz[f])
            tt = int(pt[i])
            if U <= tt or (U - tt) % 2: continue
            p = (U - tt) // 2
            if p < PMIN or len(str(p)) > NDIG: continue
            r_ = int(ridx[i]); ty, aa, bb, dd = int(typ[r_]), int(a[r_]), int(b[r_]), int(d[r_])
            if RAW:
                solutions.append((p, ty, aa, bb, dd)); continue
            nprp += 1
            if not gmpy2.is_prime(p, 8): continue
            q_, r2, s_ = p + 2 * aa, p + 2 * bb, p + 2 * dd
            if not (gmpy2.is_prime(q_, 8) and gmpy2.is_prime(r2, 8) and gmpy2.is_prime(s_, 8)): continue
            stats['four primes'] = stats.get('four primes', 0) + 1
            if gmpy2.next_prime(p) != q_ or gmpy2.next_prime(q_) != r2 or gmpy2.next_prime(r2) != s_: continue
            N = {1: p * q_ + r2 * s_, 2: p * s_ + q_ * r2, 3: p * r2 + q_ * s_}[ty]
            assert gmpy2.is_square(N)
            solutions.append((p, ty, aa, bb, dd))
            print(f"SOLUTION type {ty}: p = {p} ({len(str(p))} digits), gaps {2*aa}, {2*bb}, {2*dd}", flush=True)

import os
PART, NPART = (int(x) for x in os.environ.get('PART', '0/1').split('/'))     # run chunk c only if c % NPART == PART
t0 = time.time(); dlo = 3; chunk = 0
while dlo <= DMAX:
    dhi = dlo; est = 0
    while dhi <= DMAX and est < 1_200_000: est += dhi * dhi; dhi += 1    # about 3e6 pairs per chunk
    chunk += 1
    if chunk % NPART != PART: dlo = dhi; continue
    rec = triples(dlo, dhi)
    if rec is not None: run_chunk(rec)
    if not RAW: print(f"  d < {dhi}: {npairs} pairs, {nprp} exact tests, {len(solutions)} solutions, {time.time() - t0:.0f} s", flush=True)
    dlo = dhi
solutions = sorted(set(solutions))
print(f"pairs sieved: {npairs}; exact candidates: {nprp}; quadruples with p, q, r, s all prime: {stats.get('four primes', 0)}; solutions: {len(solutions)}")
for sol in solutions: print("   ", sol)
