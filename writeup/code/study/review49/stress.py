# Stress tests of the level (data - own law) by subsets and height bands; Hotelling and chi2 calibration
import json, numpy as np, sys
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
import analyse as A
from scipy import stats
RK = json.load(open(A.D + 'ranks.json'))
ids = sorted(i for i in A.LAW if i in A.S2 and i not in A.EXCL)
ids_all = sorted(i for i in A.LAW if i in A.S2)
pick = [1, 3, 7, 11, 15, 19]
rng = np.random.default_rng(7)


def summ(sub, band='full', law='full', lab=''):
    S = np.array([A.S2[i][band] for i in sub]); T = np.array([A.LAW[i][band][law] for i in sub]); n = len(sub)
    Dd = S - T; m = Dd.mean(0); se = Dd.std(0, ddof=1) / np.sqrt(n); z = m / se
    if n > 25:
        C = np.cov(Dd.T); T2 = n * m @ np.linalg.solve(C, m); F = T2 * (n - 20) / ((n - 1) * 20); p = stats.f.sf(F, 20, n - 20)
    else:
        p = float('nan')
    print('  %-26s n=%3d  m(l=1,2,4,6,8,10) %s  max|z| %.1f  mean z %+.2f  Hotelling p %.3g' % (lab, n, A.fmt(m[pick]), np.abs(z).max(), z.mean(), p))
    return m, se


print('level data - own law by subset (full band):')
summ(ids, lab='all 136')
summ(ids_all, lab='all 138 (dups removed)')
for r in (0, 1, 2):
    summ([i for i in ids if RK[i]['r'] == r], lab='rank %d' % r)
for w in (1, -1):
    summ([i for i in ids if RK[i]['w'] == w], lab='root number %+d' % w)
Ns = np.median([RK[i]['N'] for i in ids])
summ([i for i in ids if RK[i]['N'] < Ns], lab='N < %d' % Ns); summ([i for i in ids if RK[i]['N'] >= Ns], lab='N >= %d' % Ns)
for g in ('0', '1', '2', 'b'):
    summ([i for i in ids if A.FS[i]['g2'] == g], lab='|a2| class %s' % g)
print('\nheight bands (own law averaged over the same band):')
for b in ('lo', 'hi', 'b60'):
    summ(ids, band=b, lab='band %s' % b)
    for r in (0, 1, 2):
        summ([i for i in ids if RK[i]['r'] == r], band=b, lab='  band %s rank %d' % (b, r))
# difference lo - hi
Slo = np.array([A.S2[i]['lo'] for i in ids]) - np.array([A.LAW[i]['lo']['full'] for i in ids])
Shi = np.array([A.S2[i]['hi'] for i in ids]) - np.array([A.LAW[i]['hi']['full'] for i in ids])
d = Slo - Shi; print('  lo - hi residual: m %s, se %s' % (A.fmt(d.mean(0)[pick]), A.fmt(d.std(0, ddof=1)[pick] / np.sqrt(len(ids)), '%.4f')))
# central zero: the band 20-60 alone (law over 20-60 not computed: compare b60 with full)
dd = np.array([A.S2[i]['full'] for i in ids]) - np.array([A.LAW[i]['full']['full'] for i in ids]) - (
    np.array([A.S2[i]['b60'] for i in ids]) - np.array([A.LAW[i]['b60']['full'] for i in ids]))
print('  (full - law_full) - (b60 - law_b60): m %s se %s' % (A.fmt(dd.mean(0)[pick]), A.fmt(dd.std(0, ddof=1)[pick] / np.sqrt(len(ids)), '%.4f')))

# Hotelling calibration: bootstrap of centred residual rows
S = A.meas(ids); Tf = np.array([A.P[i]['full'] for i in ids]); T = A.arr(ids, 'full', 'full')
for lab, TT in (('frozen law', Tf), ('own law', T)):
    Dd = S - TT; n = len(Dd); m = Dd.mean(0); C = np.cov(Dd.T); T2 = n * m @ np.linalg.solve(C, m)
    E = Dd - m; bs = []
    for _ in range(4000):
        e = E[rng.integers(0, n, n)]; mm = e.mean(0); Cb = np.cov(e.T); bs.append(n * mm @ np.linalg.solve(Cb, mm))
    bs = np.array(bs); F = T2 * (n - 20) / ((n - 1) * 20)
    print('\nHotelling (%s): T2 %.1f; F-reference p %.3f; bootstrap-null p %.3f (null mean %.1f, F-ref mean %.1f)' % (
        lab, T2, stats.f.sf(F, 20, n - 20), (bs >= T2).mean(), bs.mean(), (n - 1) * 20 / (n - 20) * (n - 20) / (n - 22)))
    # eigen-decomposition: which direction carries T2
    w, V = np.linalg.eigh(C); proj = V.T @ m; contrib = n * proj ** 2 / w
    o = np.argsort(contrib)[::-1][:3]
    print('   largest T2 contributions: eigval %s, share %s' % (A.fmt(w[o], '%.2e'), A.fmt(contrib[o] / T2, '%.2f')))

# chi2 of the p = 2 / p = 3 indicator contrasts
def design(sub):
    return np.array([[1.0] + [float(A.FS[i]['g2'] == g) for g in ('0', '2', 'b')] + [float(A.FS[i]['g3'] == g) for g in ('0', '2', '3', 'b')] for i in sub])
X = design(ids); XtXi = np.linalg.pinv(X.T @ X); n, k = X.shape


def chi2s(S, TT):
    cS = XtXi @ X.T @ S; cT = XtXi @ X.T @ TT; res = S - X @ cS
    Creg = np.cov(res.T, ddof=k) * (n - 1) / (n - k)                 # as registered (residual of S itself)
    R = S - TT; cR = XtXi @ X.T @ R; resR = R - X @ cR
    Cnoise = (resR.T @ resR) / (n - k)                                # residual of S - law (the measurement noise)
    out = {}
    for lab, rows in (('p2', [1, 2, 3]), ('p3', [4, 5, 6, 7])):
        d = (cS - cT)[rows].ravel(); A_ = XtXi[np.ix_(rows, rows)]
        out[lab] = (float(d @ np.linalg.solve(np.kron(A_, Creg), d)), float(d @ np.linalg.solve(np.kron(A_, Cnoise), d)))
    return out


for lab, TT in (('frozen law', Tf), ('own law', T)):
    c0 = chi2s(S, TT)
    # null by bootstrap of the residual rows of S - law around their regression fit (H0: law correct, coefficients 0)
    R = S - TT; cR = XtXi @ X.T @ R; E = R - X @ cR
    bs = {'p2': [], 'p3': []}
    for _ in range(1000):
        Sb = TT + X[:, :1] @ cR[:1] + E[rng.integers(0, n, n)]
        cb = chi2s(Sb, TT)
        for kk in bs: bs[kk].append(cb[kk])
    for kk, cells in (('p2', 60), ('p3', 80)):
        b = np.array(bs[kk])
        print('chi2 %s (%s, %d cells): registered-cov %.1f (null mean %.1f, boot p %.3f; chi2-table p %.3g) | noise-cov %.1f (null mean %.1f, boot p %.3f; chi2-table p %.3g)' % (
            kk, lab, cells, c0[kk][0], b[:, 0].mean(), (b[:, 0] >= c0[kk][0]).mean(), stats.chi2.sf(c0[kk][0], cells),
            c0[kk][1], b[:, 1].mean(), (b[:, 1] >= c0[kk][1]).mean(), stats.chi2.sf(c0[kk][1], cells)))
