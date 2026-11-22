"""Review 26: own k-th neighbour law (power series in z, adapted from review23/kth.py), many lam at once.
Saves dist[lam, k, D] (class-averaged, per prime start) and the per-class gap table T[lam, class, D] (z^0 coefficient).
usage: mylaw.py lam1,lam2,... DMAX K out.npz"""
import sys, math, time
from itertools import combinations
import numpy as np

LAMS = [float(x) for x in sys.argv[1].split(',')]; DMAX = int(sys.argv[2]); K = int(sys.argv[3]); OUTF = sys.argv[4]
SMALL = [2, 3, 5, 7, 11, 13]; PP = 30030; DEN = math.prod(1 - 1 / q for q in SMALL)
X = 10 ** 7
s = np.ones(X + 1, bool); s[:2] = False
for i in range(2, math.isqrt(X) + 1):
    if s[i]:
        s[i * i::i] = False
PR = np.nonzero(s)[0].astype(np.float64); del s
MIDQ = [int(q) for q in PR if 13 < q <= DMAX + 1]
FARQ = PR[PR > DMAX + 1]
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
IDX = np.full(PP, -1); IDX[CLS] = np.arange(len(CLS))
DL = list(range(2, DMAX + 1, 2))
dist = np.zeros((len(LAMS), K, len(DL)))
T = np.zeros((len(LAMS), len(CLS), len(DL)), np.float32)
t0 = time.time()
for di, D in enumerate(DL):
    pts = np.arange(2, D, 2, dtype=np.int64); n = len(pts)
    S2 = S([[0, D]])[0]
    cls = CLS[np.gcd(CLS + D, PP) == 1]
    M = (np.gcd(cls[:, None] + pts[None, :], PP) == 1).astype(float)
    R = np.zeros((n, n))
    if n:
        s1 = S(np.stack([np.zeros(n, np.int64), np.full(n, D), pts], 1)) / S2
    if n > 1:
        ij = np.array(list(combinations(range(n), 2)))
        s2 = S(np.stack([np.zeros(len(ij), np.int64), np.full(len(ij), D), pts[ij[:, 0]], pts[ij[:, 1]]], 1)) / S2
        rr = s2 / (s1[ij[:, 0]] * s1[ij[:, 1]]) - 1
        R[ij[:, 0], ij[:, 1]] = rr; R[ij[:, 1], ij[:, 0]] = rr
    for li, lam in enumerate(LAMS):
        rho0 = 1 / (lam * DEN)
        poly = np.zeros((len(cls), K)); poly[:, 0] = 1
        if n:
            p = rho0 * s1
            for i in range(n):
                m = M[:, i] * p[i]
                new = poly * (1 - m)[:, None]; new[:, 1:] += poly[:, :-1] * m[:, None]; poly = new
        Qs = np.zeros((len(cls), K))
        if n > 1:
            u = p / (1 - p); a = np.zeros((K, n)); a[0] = -u
            for t in range(1, K):
                a[t] = u * (-u) ** (t - 1) * (1 + u)
            A = [M * a[t][None, :] for t in range(K)]; AR = [A[t] @ R for t in range(K)]
            for t in range(K):
                for m_ in range(t + 1):
                    Qs[:, t] += 0.5 * np.sum(AR[m_] * A[t - m_], axis=1)
        E = np.zeros_like(Qs); E[:, 0] = np.exp(Qs[:, 0])
        for t in range(1, K):
            E[:, t] = sum(mm * Qs[:, mm] * E[:, t - mm] for mm in range(1, t + 1)) / t
        G = np.zeros_like(Qs)
        for t in range(K):
            G[:, t] = sum(poly[:, m_] * E[:, t - m_] for m_ in range(t + 1))
        dist[li, :, di] = rho0 * S2 * G.sum(0) / len(CLS)
        T[li, IDX[cls], di] = rho0 * S2 * G[:, 0]
    if D % 40 == 0:
        print('D', D, round(time.time() - t0), 's', flush=True)
np.savez(OUTF, dist=dist, T=T, D=np.array(DL), lams=np.array(LAMS))
