"""Triples (0,2,6) and (0,4,6) with my estimator (review 155).  Uses my recount (T026, T046 per cell) and the regressors of
rc_fit.py (same cells).  Prediction by enumeration: zeta 3*S3, chi -> S3 * sum_i mean_{n adm mod 420} conj chi(n + b_i).
usage: rc_triple.py rc_full.npz rc_fit_out.npz
"""
import sys, math, itertools
import numpy as np
D = np.load(sys.argv[1]); F = np.load(sys.argv[2])
r, CH, NUI = F['r'], F['CH'], F['NUI']
K = np.arange(int(D['i0']) + 1, int(D['i1']) + 1); X = D['X'].astype(float); dx = X[K] - X[K - 1]; sq = np.sqrt(X[K]); u = np.log(X[K])
M = 420
# S3 for (0,2,6) (same for (0,4,6)): own product over primes up to 1e7
N = 10 ** 7
s = np.ones(N + 1, bool); s[:2] = False
for p in range(2, int(N ** 0.5) + 1):
    if s[p]: s[p * p::p] = False
P = np.nonzero(s)[0].astype(float)
P5 = P[P >= 5]
S3 = 4 * (9 / 8) * np.exp(np.sum(np.log1p(-3 / P5) - 3 * np.log1p(-1 / P5)))
print('S3 (own, primes < 1e7) = %.6f' % S3)
def dl(p, g):
    t = {}; x = 1
    for k in range(p - 1): t[x] = k; x = x * g % p
    return t
tabs = {3: dl(3, 2), 5: dl(5, 2), 7: dl(7, 3)}
def cv(p, j, c):
    if p == 4:
        if c % 2 == 0: return 0
        return 1 if j == 0 else (1 if c % 4 == 1 else -1)
    if c % p == 0: return 0
    return np.exp(2j * np.pi * j * tabs[p][c % p] / (p - 1))
chars = []
for jt in itertools.product(range(2), range(2), range(4), range(6)):
    v = np.array([cv(4, jt[0], c) * cv(3, jt[1], c) * cv(5, jt[2], c) * cv(7, jt[3], c) for c in range(M)])
    q = (4 if jt[0] else 1) * (3 if jt[1] else 1) * (5 if jt[2] else 1) * (7 if jt[3] else 1)
    chars.append((jt, v, q))
nonpr = [i for i, c in enumerate(chars) if c[2] > 1]
colinfo = []; seen = set()
for i in nonpr:
    if i in seen: continue
    jj = [k for k in nonpr if np.allclose(chars[k][1], np.conj(chars[i][1]))][0]
    seen |= {i, jj}
    colinfo += [(i, 're')] if jj == i else [(i, 're'), (i, 'im')]
def pred(bs):
    adm = [n for n in range(M) if math.gcd(n * (n + bs[1]) * (n + bs[2]), M) == 1]
    out = []
    for i, part in colinfo:
        v = chars[i][1]
        b = S3 * sum(np.mean([np.conj(v[(n + bi) % M]) for n in adm]) for bi in bs)
        out.append(b.real if part == 're' else b.imag)
    return np.array(out)
w = 1 / u ** 3
Z = np.column_stack([r, CH, NUI]); sw = np.sqrt(w); A = Z * sw[:, None]; G = np.linalg.inv(A.T @ A)
def wls(y, Z, blocks=40):
    A = Z * sw[:, None]; G = np.linalg.inv(A.T @ A); b = G @ (A.T @ (y * sw)); res = y * sw - A @ b
    s2 = res @ res / (len(y) - Z.shape[1])
    n = len(y); e = np.linspace(0, n, blocks + 1).astype(int); bs = []
    for g in range(blocks):
        m = np.ones(n, bool); m[e[g]:e[g + 1]] = False
        bs.append(np.linalg.solve(A[m].T @ A[m], A[m].T @ (y * sw)[m]))
    bs = np.array(bs); global JK
    JK = np.sqrt((blocks - 1) / blocks * ((bs - bs.mean(0)) ** 2).sum(0))
    return b, s2 * G
for name, bs in (('T026', (0, 2, 6)), ('T046', (0, 4, 6))):
    y = (D[name][K] - S3 * dx) / sq
    pv = pred(bs)
    b, cov = wls(y, Z)
    print('\n%s: zeta %.3f +- %.3f (OLS) +- %.3f (jk) (pred %.4f)' % (name, b[0], math.sqrt(cov[0, 0]), JK[0], 3 * S3))
    k3 = [k for k, (i, p) in enumerate(colinfo) if chars[i][0] == (0, 1, 0, 0)][0]
    print('  raw chi_-3 %.3f +- %.3f (pred %.4f)' % (b[1 + k3], math.sqrt(cov[1 + k3, 1 + k3]), pv[k3]))
    pres = np.abs(pv) > 1e-12
    null = ~pres
    dn = b[1:1 + CH.shape[1]][null]; Cn = cov[1:1 + CH.shape[1], 1:1 + CH.shape[1]][np.ix_(null, null)]
    print('  null set chi2 %.1f on %d' % (dn @ np.linalg.solve(Cn, dn), null.sum()))
    for label, cs in (('all present (q | 420)', None), ('q=3', {3}), ('q=5', {5}), ('q=7', {7}), ('q=15', {15}), ('q=21', {21}), ('q=35', {35}), ('q=105', {105})):
        sel = np.array([pres[k] and (cs is None or chars[colinfo[k][0]][2] in cs) for k in range(len(pv))])
        if sel.sum() == 0:
            print('  %-24s no predicted characters' % label); continue
        comb = CH[:, sel] @ pv[sel]
        rest = [k for k in range(CH.shape[1]) if not sel[k]]
        b2, c2 = wls(y, np.column_stack([r, comb, CH[:, rest], NUI]))
        print('  loudness %-24s %.3f +- %.3f (OLS) +- %.3f (jk)  (%d columns)' % (label, b2[1], math.sqrt(c2[1, 1]), JK[1], sel.sum()))
