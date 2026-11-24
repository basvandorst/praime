"""Independent recomputation (reviewer's own code) of the within part of block transition tables of consecutive primes,
plus counterfactual insertion-response counts per rough class mod 510510.
usage: mywin.py OUTDIR X0 [X0 ...]   windows [X0, X0 + 1e8)
Per window writes OUTDIR/w_<x0>.npz with
  per cell (q,k,H): per-block within/(C-1)^2 (raw, no NB/(NB-1)), per-block col/(C-1), N per block
  pooled transition counts per (q,k)
  response counts: for each rough class index u mod 510510 and q in (3,5,7): counts[u, cL, cR] for composite rough
  sites and for prime sites (cL, cR = classes of the previous and next prime)."""
import sys, os, math
import numpy as np

OUT = sys.argv[1]; WIN = [int(float(a)) for a in sys.argv[2:]]
SPAN = 10 ** 8; PAD = 4000
HS = (10 ** 4, 2 * 10 ** 4, 5 * 10 ** 4, 10 ** 5, 2 * 10 ** 5, 5 * 10 ** 5, 10 ** 6)
CELLS = ((3, 1), (3, 2), (5, 1), (7, 1), (10, 2), (10, 1))
QR = 510510

def small_primes(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.flatnonzero(s)

def sieve(lo, hi, base):
    """primes in [lo, hi), odd-only segmented sieve"""
    out = []
    SEG = 2 * 10 ** 7
    for a in range(lo, hi, SEG):
        b = min(hi, a + SEG)
        a0 = a | 1                                   # first odd >= a
        n = (b - a0 + 1) // 2                        # odd numbers a0, a0+2, ..., < b
        if n <= 0: continue
        seg = np.ones(n, bool)
        for p in base[1:].tolist():                  # skip 2
            if p * p >= b: break
            # first odd multiple of p that is >= max(a0, p*p)
            st = max(p * p, ((a0 + p - 1) // p) * p)
            if st % 2 == 0: st += p
            seg[(st - a0) // 2::p] = False
        out.append(a0 + 2 * np.flatnonzero(seg))
    P = np.concatenate(out)
    return P[P >= max(lo, 3)]

def red(q): return [a for a in range(q) if math.gcd(a, q) == 1]

def block_stats(cls, blk, NB, C, k):
    """cls: class index per prime (-1 if not reduced), blk: block of p_n.  Returns per-block within, col, N."""
    a = cls[:-k]; b = cls[k:]; bl = blk[:-k]
    ok = (bl >= 0) & (bl < NB) & (a >= 0) & (b >= 0)
    idx = (bl[ok] * C + a[ok]) * C + b[ok]
    X = np.bincount(idx, minlength=NB * C * C).reshape(NB, C, C).astype(float)
    Na = X.sum(2)
    T = X.sum(0); Ph = T / T.sum(1, keepdims=True)
    E = Na[:, :, None] * Ph[None]
    chi = ((X - E) ** 2 / E).sum((1, 2))
    cv = (X - E).sum(1)                               # column sums, sum to 0
    # covariance of column sums under row-multinomial: M = sum_a Na (diag Ph_a - Ph_a Ph_a^T); drop the last coordinate
    Mfull = np.einsum('ba,aij->bij', Na, np.stack([np.diag(Ph[i]) - np.outer(Ph[i], Ph[i]) for i in range(C)]))
    Mr = Mfull[:, :-1, :-1]; cr = cv[:, :-1]
    sol = np.linalg.solve(Mr, cr[:, :, None])[:, :, 0]
    col = (cr * sol).sum(1)
    return chi - col, col, Na.sum(1), T

base = small_primes(int(math.isqrt(max(WIN) + SPAN + PAD)) + 2)
# rough classes mod QR
ok = np.ones(QR, bool)
for p in (2, 3, 5, 7, 11, 13, 17): ok[::p] = False
RIDX = -np.ones(QR, np.int64); RU = np.flatnonzero(ok); RIDX[RU] = np.arange(len(RU)); NR = len(RU)

for x0 in WIN:
    P = sieve(x0 - PAD, x0 + SPAN + PAD, base)
    res = {}
    for q, k in CELLS:
        r = red(q); C = len(r); lut = -np.ones(q, np.int64); lut[r] = np.arange(C); cls = lut[P % q]
        for H in HS:
            NB = SPAN // H
            blk = np.where(P >= x0, (P - x0) // H, -1)
            w, c, N, T = block_stats(cls, blk, NB, C, k)
            res['win_%d_%d_%d' % (q, k, H)] = w / (C - 1) ** 2
            res['col_%d_%d_%d' % (q, k, H)] = c / (C - 1)
            if H == HS[0]:
                res['N_%d_%d' % (q, k)] = N; res['T_%d_%d' % (q, k)] = T
    # insertion-response counts per rough class mod QR (k = 1), sites in [x0, x0 + SPAN)
    for q in (3, 5, 7):
        r = red(q); C = len(r); lut = -np.ones(q, np.int64); lut[r] = np.arange(C)
        res['rc_comp_%d' % q] = np.zeros(NR * C * C, np.int64); res['rc_prime_%d' % q] = np.zeros(NR * C * C, np.int64)
    CH = 10 ** 7
    for a in range(x0, x0 + SPAN, CH):
        xs = a + np.arange(CH, dtype=np.int64)
        xs = xs[ok[xs % QR]]
        j = np.searchsorted(P, xs)                   # P[j] >= x
        isp = P[j] == xs
        prevp = P[j - 1]; nextp = np.where(isp, P[j + 1], P[j])
        u = RIDX[xs % QR]
        for q in (3, 5, 7):
            r = red(q); C = len(r); lut = -np.ones(q, np.int64); lut[r] = np.arange(C)
            cL = lut[prevp % q]; cR = lut[nextp % q]
            idx = (u * C + cL) * C + cR
            res['rc_comp_%d' % q] += np.bincount(idx[~isp], minlength=NR * C * C)
            res['rc_prime_%d' % q] += np.bincount(idx[isp], minlength=NR * C * C)
    res['nprimes'] = np.array([np.sum((P >= x0) & (P < x0 + SPAN))]); res['x0'] = np.array([x0])
    res['first_last'] = np.array([P[P >= x0][0], P[P < x0 + SPAN][-1]])
    np.savez_compressed(os.path.join(OUT, 'w_%d.npz' % x0), **res)
    print('window %d done: %d primes' % (x0, res['nprimes'][0]), flush=True)
