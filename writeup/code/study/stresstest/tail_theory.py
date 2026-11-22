"""Pre-registration 3: long-gap cells (max(g1,g2) in 38..70), G = 140, both heights. Same law construction as theory.py."""
import sys, math, datetime
sys.path.insert(0, '/home/bas/praime/writeup/code/study')
import numpy as np, r47
PP = r47.PP
CLS = np.array([a for a in range(PP) if math.gcd(a, PP) == 1], dtype=np.int64)
IDX = np.full(PP, -1); IDX[CLS] = np.arange(len(CLS)); NC = len(CLS)
G = 140; NG = G // 2
BINS = [(2, 36), (38, 52), (54, 70)]
def run(lam):
    RHO0 = 1 / (lam * r47.RHO_DEN)
    P1 = np.zeros((NC, NG)); T = np.zeros((NG, NG))
    for g in range(2, G + 1, 2):
        st = r47.setup(g); pts = np.arange(2, g, 2, dtype=np.int64); n = len(pts)
        rows = IDX[CLS[np.gcd(CLS + g, PP) == 1]]
        M = st['M']; p = RHO0 * st['s1']
        L = M @ np.log1p(-p) if n else np.zeros(len(rows))
        if n > 1:
            w = p / (1 - p); V = M @ ((w[:, None] * w[None, :]) * st['R']); Qd = 0.5 * np.sum(V * M, axis=1)
        else:
            V = np.zeros((len(rows), n)); Qd = np.zeros(len(rows))
        P1[rows, g // 2 - 1] = RHO0 * st['S2'] * np.exp(L + Qd)
        for k in range(n):
            g1 = int(pts[k])
            if g1 > 70 or g - g1 > 70: continue
            ok = M[:, k] == 1
            S3 = float(r47.S_gt(np.array([[0, g1, g]]))[0]); oth = np.delete(np.arange(n), k)
            if len(oth):
                s1p = r47.S_gt(np.stack([np.zeros(len(oth), np.int64), np.full(len(oth), g1), np.full(len(oth), g), pts[oth]], 1)) / S3
                Lp = M[ok][:, oth] @ np.log1p(-RHO0 * s1p)
            else:
                Lp = np.zeros(int(ok.sum()))
            T[g1 // 2 - 1, (g - g1) // 2 - 1] = (RHO0 ** 2 * S3 * np.exp(Lp + Qd[ok] - V[ok, k])).sum() / NC
    Tc = np.zeros((NG, NG))
    for j in range(NG):
        g1 = 2 * (j + 1); nxt = IDX[(CLS + g1) % PP]; okc = nxt >= 0
        Tc[j, :] = (P1[okc, j][:, None] * P1[nxt[okc], :]).sum(0) / NC
    return T[:35, :35], Tc[:35, :35]
def cells(A):
    C = np.zeros((3, 3)); C6 = np.zeros((3, 3, 2))
    for a in range(2, 71, 2):
        for b in range(2, 71, 2):
            i = [k for k, (lo, hi) in enumerate(BINS) if lo <= a <= hi][0]; j = [k for k, (lo, hi) in enumerate(BINS) if lo <= b <= hi][0]
            C[i, j] += A[a // 2 - 1, b // 2 - 1]; C6[i, j, int(a % 6 == 0 and b % 6 == 0)] += A[a // 2 - 1, b // 2 - 1]
    return C, C6
out = {}
f = open('/home/bas/scratch/stress-B/prediction_tail.txt', 'w')
f.write('PREDICTION written %s UTC before counting long-gap cells. Cells: g1 bin x g2 bin, bins %s, excluding (A,A); plus each cell split by (6|g1 and 6|g2). Probabilities x 1e4.\n' % (datetime.datetime.now(datetime.timezone.utc).isoformat(), BINS))
for tag, lam in (('1e9', math.log(1.05e9)), ('1e12', math.log(1e12 + 5e7))):
    T, Tc = run(lam); out[tag + '_T'] = T; out[tag + '_Tc'] = Tc
    for nm, A in (('LAW', T), ('CHAIN', Tc)):
        C, C6 = cells(A)
        f.write('%s %s cells:\n%s\n both-mult-6 part:\n%s\n' % (tag, nm, np.array2string(1e4 * C, precision=4), np.array2string(1e4 * C6[:, :, 1], precision=4)))
f.close()
np.savez('/home/bas/scratch/stress-B/tail.npz', **out)
print(open('/home/bas/scratch/stress-B/prediction_tail.txt').read())
