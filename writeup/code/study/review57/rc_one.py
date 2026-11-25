"""Print the coefficients of two named odd characters of conductor 105 (the largest null deviations on the holdout) on a range."""
import sys, math, itertools
import numpy as np
D = np.load(sys.argv[1]); F = np.load(sys.argv[2])
r, CH, NUI, w = F['r'], F['CH'], F['NUI'], F['w']
K = np.arange(int(D['i0']) + 1, int(D['i1']) + 1); X = D['X'].astype(float); dx = X[K] - X[K - 1]; sq = np.sqrt(X[K])
M = 420; S = 2 * 0.6601618158468695739
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
    colinfo += [(chars[i][0], 're')] if jj == i else [(chars[i][0], 're'), (chars[i][0], 'im')]
Z = np.column_stack([r, CH, NUI]); sw = np.sqrt(w); A = Z * sw[:, None]; G = np.linalg.inv(A.T @ A)
for h in (2, 4):
    y = (D['P%d' % h][K] - S * dx) / sq
    b = G @ (A.T @ (y * sw)); res = y * sw - A @ b; s2 = res @ res / (len(y) - Z.shape[1]); se = np.sqrt(s2 * np.diag(G))
    for target in [((0, 1, 2, 2), 'im'), ((0, 1, 1, 5), 'im'), ((0, 0, 1, 0), 're'), ((0, 0, 1, 0), 'im')]:
        k = colinfo.index(target)
        print('h=%d %s %s: %+.4f +- %.4f (z %+.2f)' % (h, target[0], target[1], b[1 + k], se[1 + k], b[1 + k] / se[1 + k]))
