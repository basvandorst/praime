"""Independent k-th neighbour law (nothing imported from the study). Own S_{>13}; generating function expanded as a
power series in z (no Fourier sum): G = prod_i (1 - p_i + p_i z) * exp(Q(z)), Q = 1/2 sum_{i!=j} R_ij w_i w_j,
w_i = p_i (z-1)/(1 + p_i (z-1)) = sum_n a_{i,n} z^n.
usage: kth.py lam DMAX K q1,q2,.. out.npz [Dlist]"""
import sys, math, time
from itertools import combinations
import numpy as np

LAM = float(sys.argv[1]); DMAX = int(sys.argv[2]); K = int(sys.argv[3]); QS = [int(x) for x in sys.argv[4].split(',')]
OUTF = sys.argv[5]
DL = [int(x) for x in sys.argv[6].split(',')] if len(sys.argv) > 6 else list(range(2, DMAX + 1, 2))
SMALL = [2, 3, 5, 7, 11, 13]; PP = 30030
RHO0 = 1 / (LAM * math.prod(1 - 1 / q for q in SMALL))
X = 10 ** 7
s = np.ones(X + 1, bool); s[:2] = False
for i in range(2, math.isqrt(X) + 1):
    if s[i]:
        s[i * i::i] = False
PR = np.nonzero(s)[0].astype(np.float64)
MIDQ = [int(q) for q in PR if 13 < q <= DMAX + 1]
FARQ = PR[PR > DMAX + 1]
# for q > DMAX+1 all points of a set (diameter <= DMAX) are distinct mod q: factor (1 - k/q)/(1 - 1/q)^k
FAR = {k: math.exp(float(np.sum(np.log1p(-k / FARQ) - k * np.log1p(-1 / FARQ))) - math.comb(k, 2) / (X * math.log(X)))
       for k in range(1, 5)}


def S(sets):
    sets = np.asarray(sets, np.int64); N, k = sets.shape
    lg = np.full(N, math.log(FAR[k]))
    for q in MIDQ:
        r = np.sort(sets % q, axis=1); nu = 1 + np.count_nonzero(np.diff(r, axis=1), axis=1)
        lg += np.log(1 - nu / q) - k * math.log(1 - 1 / q)
    return np.exp(lg)


CLS = np.array([a for a in range(PP) if math.gcd(a, PP) == 1], np.int64)
dist = np.zeros((K, len(DL)))
gap_per_class = {}
t0 = time.time()
for di, D in enumerate(DL):
    pts = np.arange(2, D, 2, dtype=np.int64); n = len(pts)
    S2 = S([[0, D]])[0]
    cls = CLS[np.gcd(CLS + D, PP) == 1]
    M = (np.gcd(cls[:, None] + pts[None, :], PP) == 1).astype(float)
    if n:
        s1 = S(np.stack([np.zeros(n, np.int64), np.full(n, D), pts], 1)) / S2
        p = RHO0 * s1
    R = np.zeros((n, n))
    if n > 1:
        ij = np.array(list(combinations(range(n), 2)))
        s2 = S(np.stack([np.zeros(len(ij), np.int64), np.full(len(ij), D), pts[ij[:, 0]], pts[ij[:, 1]]], 1)) / S2
        rr = s2 / (s1[ij[:, 0]] * s1[ij[:, 1]]) - 1
        R[ij[:, 0], ij[:, 1]] = rr; R[ij[:, 1], ij[:, 0]] = rr
    # polynomial part, coefficients 0..K-1
    poly = np.zeros((len(cls), K)); poly[:, 0] = 1
    for i in range(n):
        m = M[:, i] * p[i]
        new = poly * (1 - m)[:, None]
        new[:, 1:] += poly[:, :-1] * m[:, None]
        poly = new
    # quadratic part as a series
    Qs = np.zeros((len(cls), K))
    if n > 1:
        u = p / (1 - p)
        a = np.zeros((K, n)); a[0] = -u
        for t in range(1, K):
            a[t] = u * (-u) ** (t - 1) * (1 + u)
        A = [M * a[t][None, :] for t in range(K)]
        AR = [A[t] @ R for t in range(K)]
        for t in range(K):
            for m in range(t + 1):
                Qs[:, t] += 0.5 * np.sum(AR[m] * A[t - m], axis=1)
    E = np.zeros_like(Qs); E[:, 0] = np.exp(Qs[:, 0])
    for t in range(1, K):
        E[:, t] = sum(mm * Qs[:, mm] * E[:, t - mm] for mm in range(1, t + 1)) / t
    G = np.zeros_like(Qs)
    for t in range(K):
        G[:, t] = sum(poly[:, m] * E[:, t - m] for m in range(t + 1))
    dist[:, di] = RHO0 * S2 * G.sum(0) / len(CLS)
    gap_per_class[D] = (cls, RHO0 * S2 * G[:, 0])
    if D % 60 == 0:
        print('D', D, round(time.time() - t0), 's', flush=True)
np.savez(OUTF, dist=dist, D=np.array(DL), lam=LAM)
if len(DL) == DMAX // 2:
    Ds = np.array(DL)
    for k in range(1, K + 1):
        dk = dist[k - 1]
        print(k, 'mass %.6f' % dk.sum(), ' '.join('%d:%.6f' % (q, dk[Ds % q == 0].sum() / dk.sum()) for q in QS), flush=True)
    # per-class gap law table for the class-chain comparator of round 83
    T = np.zeros((len(CLS), len(DL))); IDX = np.full(PP, -1); IDX[CLS] = np.arange(len(CLS))
    for di, D in enumerate(DL):
        c, v = gap_per_class[D]; T[IDX[c], di] = v
    np.save(OUTF.replace('.npz', '_gaptable.npy'), T)
