"""Derived index sequences of the Mersenne exponents (real and null), r4_index experiment.

Real: the 51 exponents below B = 83,880,077 (GIMPS: complete below B).
Primes below B: segmented sieve of the odd numbers (segments of 10^7).

Derived sequences (each a list of positive integers):
  pi    : pi(a_n)                                          n = 1..51
  gap   : pi(a_{n+1}) - pi(a_n)                            50 terms
  idx1  : index of a_n among primes = 1 mod 4, for the a_n = 1 mod 4
  idx3  : index of a_n among primes = 3 mod 4, for the a_n = 3 mod 4
  cls   : index of a_n among primes in its own class mod 4 (n >= 2; interleaved)
  tf    : index among primes p < XTF = 10^6 with no factor 2jp+1 of 2^p-1 (j <= 1000, factor < 2^p-1)
  blen  : a_n - 2^(L-1) + 1, L = bit length (index among integers of the same bit length)
  pblen : index among primes of the same bit length
  pdec  : index among primes with the same number of decimal digits
In-sample flag of a term: its exponent is among a_1..a_25 (gap n: both a_n, a_{n+1} among a_1..a_25).

Nulls:
  W  (all sequences but tf): each prime p < B taken with probability min(1, e^gamma log2(a p)/p),
     a = 6 for p = 1 mod 4, a = 2 otherwise; conditioned on exactly 51 terms (rejection).
     Exact sampler: p < 10^5 direct Bernoulli; above, blocks of 4096 primes, Bernoulli(qmax_block) by a
     binomial count + uniform distinct positions, then thinning by q/qmax.
  TF (tf only): each survivor p < 10^6 taken with probability min(1, e^gamma log2(2000 p)/p)
     (heuristic P(M_p prime | no factor below Q) = e^gamma log2(Q)/p with Q = 2000p, MERSENNE.md sec. 2),
     conditioned on the real count of exponents below 10^6 (33).
Output: data.npz
"""
import numpy as np, math, sys, time

B = 83_880_077
XTF, JTF = 10 ** 6, 1000
REAL = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253, 4423,
        9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243, 110503, 132049, 216091, 756839, 859433,
        1257787, 1398269, 2976221, 3021377, 6972593, 13466917, 20996011, 24036583, 25964951, 30402457,
        32582657, 37156667, 42643801, 43112609, 57885161, 74207281, 77232917, 82589933]
assert len(REAL) == 51 and REAL == sorted(REAL) and REAL[-1] < B
EG = math.exp(0.5772156649015329)
NIN = 25
NNULL = int(sys.argv[1]) if len(sys.argv) > 1 else 300


def seg_sieve(n, seg=10 ** 7):
    r = int(math.isqrt(n)) + 1
    s = np.ones(r + 1, bool); s[:2] = False
    for p in range(2, int(math.isqrt(r)) + 1):
        if s[p]: s[p * p::p] = False
    base = np.nonzero(s)[0][1:]  # odd base primes
    out = [np.array([2], np.uint32)]
    for lo in range(1, n, seg):   # odd numbers in [lo, hi)
        hi = min(lo + seg, n)
        m = np.ones((hi - lo + 1) // 2, bool)  # m[i] <-> lo + 2i
        for p in base:
            if p * p >= hi: break
            st = max(p * p, ((lo + p - 1) // p) * p)
            if st % 2 == 0: st += p
            m[(st - lo) // 2::p] = False
        v = lo + 2 * np.nonzero(m)[0]
        v = v[(v >= 3) & (v < n)]
        out.append(v.astype(np.uint32))
    return np.concatenate(out)


t0 = time.time()
PR = seg_sieve(B)
print('primes below B:', len(PR), 'time %.1f' % (time.time() - t0))
assert len(PR) == 4_882_467
c1 = np.cumsum(PR % 4 == 1).astype(np.int64)
c3 = np.cumsum(PR % 4 == 3).astype(np.int64)
P64 = PR.astype(np.int64)

# ---- trial-factoring survivors below XTF (q = 2jp+1, j <= JTF, q = +-1 mod 8, proper factor) ----
t0 = time.time()
P = P64[(P64 >= 3) & (P64 < XTF)].astype(np.uint64)
assert 2 * JTF * int(P.max()) + 1 < 2 ** 31
kmin = np.zeros(len(P), np.int64)
nb = int(P.max()).bit_length()
for k in range(1, JTF + 1):
    open_ = kmin == 0
    p = P[open_]; q = 2 * np.uint64(k) * p + np.uint64(1)
    ok = (q % np.uint64(8) == 1) | (q % np.uint64(8) == 7)
    p, q = p[ok], q[ok]
    r = np.ones(len(q), np.uint64)
    for b in range(nb - 1, -1, -1):
        r = (r * r) % q
        bit = (p >> np.uint64(b)) & np.uint64(1)
        r = np.where(bit == 1, (r * np.uint64(2)) % q, r)
    hit = r == 1
    small = p < 64
    if small.any():
        proper = np.array([int(qq) < (1 << int(pp)) - 1 if pp < 64 else True for pp, qq in zip(p, q)])
        hit &= proper
    kmin[np.flatnonzero(open_)[ok][hit]] = k
SURV = np.concatenate([[2], P[kmin == 0].astype(np.int64)])
print('TF survivors below %d (j <= %d): %d of %d primes, time %.1f' % (XTF, JTF, len(SURV), len(P) + 1, time.time() - t0))
assert all(e in set(SURV.tolist()) for e in REAL if e < XTF)


def derive(a):
    """a: sorted list of primes < B. Returns dict name -> (values list, in-sample flags list)."""
    a = np.asarray(a, np.int64)
    n = len(a)
    ins = np.arange(n) < NIN
    pos = np.searchsorted(P64, a)
    assert (P64[pos] == a).all()
    pi = pos + 1
    D = {}
    D['pi'] = (pi, ins)
    D['gap'] = (np.diff(pi), ins[1:])
    m1 = a % 4 == 1; m3 = a % 4 == 3
    D['idx1'] = (c1[pos][m1], ins[m1])
    D['idx3'] = (c3[pos][m3], ins[m3])
    cl = np.where(m1, c1[pos], c3[pos])
    D['cls'] = (cl[a > 2], ins[a > 2])
    L = np.array([int(x).bit_length() for x in a])
    D['blen'] = (a - (1 << (L - 1)) + 1, ins)
    D['pblen'] = (pi - np.searchsorted(P64, 1 << (L - 1), 'left'), ins)
    dd = np.array([len(str(int(x))) for x in a])
    D['pdec'] = (pi - np.searchsorted(P64, 10 ** (dd - 1), 'left'), ins)
    return D


def derive_tf(a):
    a = np.asarray(a, np.int64)
    a = a[a < XTF]
    pos = np.searchsorted(SURV, a)
    assert (SURV[pos] == a).all()
    return pos + 1, np.arange(len(a)) < NIN


NAMES = ['pi', 'gap', 'idx1', 'idx3', 'cls', 'tf', 'blen', 'pblen', 'pdec']

# ---- null W ----
rng = np.random.default_rng(20261115)
pf = P64.astype(np.float64)
aa = np.where(P64 % 4 == 1, 6.0, 2.0)
QW = np.minimum(1.0, EG * np.log2(aa * pf) / pf)
print('expected count below B under W: %.2f' % QW.sum())
K0 = int(np.searchsorted(P64, 10 ** 5))
BS = 4096
starts = np.arange(K0, len(P64), BS)
ends = np.minimum(starts + BS, len(P64))
qmax = np.maximum.reduceat(QW[K0:], starts - K0)
nblk = ends - starts


def draw_W():
    sel = [np.nonzero(rng.random(K0) < QW[:K0])[0]]
    cnt = rng.binomial(nblk, qmax)
    for bi in np.nonzero(cnt)[0]:
        while True:
            posb = rng.integers(0, nblk[bi], cnt[bi])
            if len(np.unique(posb)) == cnt[bi]: break
        idx = starts[bi] + posb
        keep = rng.random(len(idx)) < QW[idx] / qmax[bi]
        sel.append(idx[keep])
    idx = np.sort(np.concatenate(sel))
    return P64[idx]


def draw_cond(drawer, n, want):
    out, tries = [], 0
    while len(out) < n:
        s = drawer(); tries += 1
        if len(s) == want:
            out.append(s)
    return out, tries


t0 = time.time()
NW, triesW = draw_cond(draw_W, NNULL, len(REAL))
print('null W: %d lists of 51 from %d draws (%.1f s)' % (len(NW), triesW, time.time() - t0))

QT = np.minimum(1.0, EG * np.log2(2000.0 * SURV) / SURV)
nreal_tf = sum(e < XTF for e in REAL)
print('null TF: expected count below %d = %.2f, real %d' % (XTF, QT.sum(), nreal_tf))


def draw_T():
    return SURV[rng.random(len(SURV)) < QT]


t0 = time.time()
NT, triesT = draw_cond(draw_T, NNULL, nreal_tf)
print('null TF: %d lists of %d from %d draws (%.1f s)' % (len(NT), nreal_tf, triesT, time.time() - t0))

# ---- assemble padded matrices: row 0 = real, rows 1.. = nulls ----
out = {}
for nm in NAMES:
    rows = []
    if nm == 'tf':
        rows.append(derive_tf(REAL))
        rows += [derive_tf(s) for s in NT]
    else:
        rows.append(derive(REAL)[nm])
        rows += [derive(s)[nm] for s in NW]
    Lmax = max(len(v) for v, _ in rows)
    M = np.full((len(rows), Lmax), -2, np.int64)
    I = np.zeros((len(rows), Lmax), bool)
    ln = np.zeros(len(rows), np.int64)
    for i, (v, f) in enumerate(rows):
        M[i, :len(v)] = v; I[i, :len(v)] = f; ln[i] = len(v)
    assert (M[M != -2] >= 1).all()
    out[nm] = M; out[nm + '_ins'] = I; out[nm + '_len'] = ln
    print(f'{nm:6s} real ({ln[0]} terms, {I[0].sum()} in-sample): {list(M[0, :ln[0]])}')
    print(f'        null lengths median {np.median(ln[1:]):.0f} range {ln[1:].min()}-{ln[1:].max()}; example {list(M[1, :12])}')
np.savez('data.npz', names=np.array(NAMES), real=np.array(REAL), surv=SURV,
         nullW=np.array(NW), nullT=np.array([np.pad(s, (0, 0)) for s in NT]), **out)
print('saved data.npz')
