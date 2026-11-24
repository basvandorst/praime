"""Internal check of the class-response decomposition (derivation.md section 6.4), no primes involved.
Model R_Q: the sieve model M_Q (Q = 30030) with every class count mod Q FIXED in every block of 33 periods
(990,990 integers): in each block, class c gets exactly K_c 'primes' at random among its 33 sites, K_c = floor(33p)
or +1 by a fixed per-class pattern. All block-level class-count fluctuations are then zero, so the law predicts
  E chi^2(R_Q) = E chi^2(M_Q) - FO_M        (FO_HL replaced by 0).
This script measures E chi^2(R_Q) for the block-homogeneity statistics; simlaw.py (same lam, same block span) gives
E chi^2(M_Q) and FO_M. usage: rigid_check.py x nchunks blocks_per_chunk seed out.json"""
import sys, math, json, time
import numpy as np

X = float(sys.argv[1]); NCH = int(sys.argv[2]); BPC = int(sys.argv[3]); SEED = int(sys.argv[4]); OUT = sys.argv[5]
LAM = math.log(X); Q = 30030
adm = np.array([a for a in range(Q) if math.gcd(a, Q) == 1], np.int64); NC = len(adm)
p = (Q / NC) / LAM; PER = 33
rng = np.random.default_rng(SEED)
K = np.floor(PER * p).astype(int) + (rng.random(NC) < (PER * p - math.floor(PER * p)))
SPECS = [('trans', 3, 1), ('trans', 3, 2), ('trans', 10, 1), ('trans', 10, 2), ('band', 0, 4)]


def cell(spec, left, right):
    kind, q, k = spec
    d = right - left
    if kind == 'band':
        c = (d - 2) // 12
        return np.where((d <= 240) & (c >= 0) & (c < 20), c, -1)
    red = [a for a in range(q) if math.gcd(a, q) == 1]
    lut = np.full(q, -1); lut[red] = np.arange(len(red)); R = len(red)
    return lut[left % q] * R + lut[right % q]


counts = {s: [] for s in SPECS}
t0 = time.time()
for ch in range(NCH):
    keys = rng.random((BPC, PER, NC))
    rank = np.argsort(np.argsort(keys, axis=1), axis=1)
    isp = rank < K[None, None, :]
    base = (np.arange(BPC * PER, dtype=np.int64)[:, None] * Q + adm[None, :]).reshape(BPC, PER, NC) + Q
    P = np.sort(base[isp]); del keys, rank, base
    blk = (P - Q) // (PER * Q)
    for s in SPECS:
        k = s[2]
        c = cell(s, P[:-k], P[k:]); b = blk[:-k]
        ncell = 20 if s[0] == 'band' else int((len([a for a in range(s[1]) if math.gcd(a, s[1]) == 1])) ** 2)
        sel = c >= 0
        counts[s].append(np.bincount(b[sel] * ncell + c[sel], minlength=BPC * ncell).reshape(BPC, ncell))
    print('chunk', ch, '%.0f s' % (time.time() - t0), flush=True)
res = {'x': X, 'lam': LAM, 'p': p, 'blocks': NCH * BPC, 'stats': {}}
for s in SPECS:
    C = np.concatenate(counts[s]).astype(float); NB = C.shape[0]; pooled = C.sum(0)
    if s[0] == 'band':
        Pc = pooled / pooled.sum(); E = C.sum(1, keepdims=True) * Pc[None, :]; dof = 19
    else:
        R = int(round(math.sqrt(C.shape[1]))); rows = np.arange(C.shape[1]) // R
        rowtot = np.zeros((NB, R)); np.add.at(rowtot.T, rows, C.T)
        prow = np.bincount(rows, weights=pooled, minlength=R); Pc = pooled / prow[rows]
        E = rowtot[:, rows] * Pc[None, :]; dof = R * (R - 1)
    chi = ((C - E) ** 2 / E).sum(1)
    res['stats']['%s:%d:%d' % s] = {'chi2_rigid': float(chi.mean() * NB / (NB - 1)), 'err': float(chi.std(ddof=1) / math.sqrt(NB)), 'dof': dof}
    print(s, res['stats']['%s:%d:%d' % s], flush=True)
json.dump(res, open(OUT, 'w'), indent=1)
