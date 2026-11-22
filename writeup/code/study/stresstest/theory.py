"""Stress test B, part 1 (theory only, no prime counting): the two-gap law of r79.py, kept per (g1, g) instead of
summed by span, for g = g1 + g2 <= GMAX. Also the chain of one-gap laws and independent gaps (law marginals).
Writes theory.npz and prediction.txt (with a timestamp) BEFORE any data are counted."""
import sys, os, math, time, datetime, json
sys.path.insert(0, '/home/bas/praime/writeup/code/study')
import numpy as np
import r47

PP = r47.PP
CLS = np.array([a for a in range(PP) if math.gcd(a, PP) == 1], dtype=np.int64)
IDX = np.full(PP, -1); IDX[CLS] = np.arange(len(CLS))
NC = len(CLS)
LAM = math.log(1.05e9); RHO0 = 1 / (LAM * r47.RHO_DEN)
G = 76; NG = G // 2          # spans g <= 76 suffice for g1, g2 <= 36 (and the mod-6 table restricted the same way)

t0 = time.time()
P1 = np.zeros((NC, NG))                  # P_a(gap = g)
T = np.zeros((NG, NG))                   # T[g1, g2] = mean_a P_a(next at g1, 2nd at g1+g2)  (index g/2 - 1)
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
        g1 = int(pts[k]); ok = M[:, k] == 1
        S3 = float(r47.S_gt(np.array([[0, g1, g]]))[0])
        oth = np.delete(np.arange(n), k)
        if len(oth):
            s1p = r47.S_gt(np.stack([np.zeros(len(oth), np.int64), np.full(len(oth), g1), np.full(len(oth), g), pts[oth]], 1)) / S3
            Lp = M[ok][:, oth] @ np.log1p(-RHO0 * s1p)
        else:
            Lp = np.zeros(int(ok.sum()))
        val = RHO0 ** 2 * S3 * np.exp(Lp + Qd[ok] - V[ok, k])
        T[g1 // 2 - 1, (g - g1) // 2 - 1] = val.sum() / NC
print('law done', round(time.time() - t0), 's', flush=True)

Tc = np.zeros((NG, NG))                  # chain: mean_a P_a(g1) P_{a+g1}(g2)
for j in range(NG):
    g1 = 2 * (j + 1); nxt = IDX[(CLS + g1) % PP]; okc = nxt >= 0
    Tc[j, :] = (P1[okc, j][:, None] * P1[nxt[okc], :]).sum(0) / NC
m1 = P1.mean(0)                          # law marginal P(gap = g)
Ti = np.outer(m1, m1)

# cells, fixed before data: g1, g2 in 2..36, size bins
BINS = [(2, 4), (6, 8), (10, 12), (14, 18), (20, 26), (28, 36)]
def binof(g):
    for i, (lo, hi) in enumerate(BINS):
        if lo <= g <= hi:
            return i
    return -1
gs = np.arange(2, 37, 2)
def cells(A):
    C = np.zeros((len(BINS), len(BINS)))
    for g1 in gs:
        for g2 in gs:
            C[binof(g1), binof(g2)] += A[g1 // 2 - 1, g2 // 2 - 1]
    return C
def mod6(A):
    C = np.zeros((3, 3))
    for g1 in gs:
        for g2 in gs:
            C[(g1 % 6) // 2, (g2 % 6) // 2] += A[g1 // 2 - 1, g2 // 2 - 1]
    return C
np.savez('/home/bas/scratch/stress-B/theory.npz', T=T, Tc=Tc, Ti=Ti, P1m=m1)
stamp = datetime.datetime.now(datetime.timezone.utc).isoformat()
with open('/home/bas/scratch/stress-B/prediction.txt', 'w') as f:
    f.write('PREDICTION written %s UTC, before any primes in [1e9, 1.1e9) were counted.\n' % stamp)
    f.write('Quantity: P(p_{n+1}-p_n = g1, p_{n+2}-p_{n+1} = g2 | p_n prime in [1e9,1.1e9)), lambda = ln(1.05e9), y = 13 conditioned law (r79 construction).\n')
    f.write('Expected counts = N_primes_with_two_successors * probability. Test: Pearson chi2 over the 36 size cells (df 36, no fitted parameter) and the 9 mod-6 cells (g1,g2<=36).\n')
    f.write('Law total mass g1,g2<=36: %.6f   chain: %.6f   indep: %.6f\n' % (cells(T).sum(), cells(Tc).sum(), cells(Ti).sum()))
    for name, A in (('LAW', T), ('CHAIN', Tc), ('INDEP', Ti)):
        f.write('\n%s size cells (rows g1 bin, cols g2 bin; bins %s), probability x 1e4:\n' % (name, BINS))
        for r in cells(A):
            f.write(' '.join('%9.4f' % (1e4 * x) for x in r) + '\n')
        f.write('%s mod-6 cells (rows g1 mod 6 = 0,2,4; cols g2 mod 6), probability x 1e4:\n' % name)
        for r in mod6(A):
            f.write(' '.join('%9.4f' % (1e4 * x) for x in r) + '\n')
    # (b)-type conditional and a few single cells
    for name, A in (('LAW', T), ('CHAIN', Tc), ('INDEP', Ti)):
        M6 = mod6(A)
        f.write('%s P(g2 = 0 mod 6 | g1 = 0 mod 6) = %.5f ; P(g2=0 mod 6 | g1 != 0 mod 6) = %.5f\n' % (name, M6[0, 0] / M6[0].sum(), M6[1:, 0].sum() / M6[1:].sum()))
    for g1, g2 in ((2, 4), (4, 2), (2, 10), (6, 6), (12, 12), (30, 6), (6, 30), (2, 2), (4, 8)):
        f.write('cell (%d,%d): law %.6e chain %.6e indep %.6e\n' % (g1, g2, T[g1 // 2 - 1, g2 // 2 - 1], Tc[g1 // 2 - 1, g2 // 2 - 1], Ti[g1 // 2 - 1, g2 // 2 - 1]))
print(open('/home/bas/scratch/stress-B/prediction.txt').read())
