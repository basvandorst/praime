"""Where does the variance contrast live? y = log|L'| - ln L = near_k + far_k, near_k = sum over the k nearest zeros on
each side of ln|unfolded distance|. Contrast per prime of Var(near), Var(far), 2 Cov(near, far) (registered design).
Also: the remove-the-wave map applied to near and far separately (does the beyond-closure part sit in near or far?)."""
import os, numpy as np, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from myfit import load, design, ols, PR, SETS
R = os.path.dirname(os.path.abspath(__file__))+'/../../../../'
rows = [r for r in load() if r['D'] not in (197, -395)]
Z = {s: np.load(fz) for s, (fv, fz) in SETS.items()}
X, nn = design(rows)
res = {}
for k in (1, 2):
    vn, vf, cv, vy = [], [], [], []
    for r in rows:
        z = np.sort(Z[r['set']]['D%d' % r['D']]); z = z[z > 0]
        idx = np.searchsorted(z, r['g'] - 1e-9)
        ok = (idx - k >= 0) & (idx + k < len(z))
        g = r['g'][ok]; L = r['Lj'][ok]; i0 = idx[ok]
        y = r['la'][ok] - np.log(L)
        near = np.zeros(len(g))
        for j in range(1, k + 1):
            near += np.log(np.abs(z[i0 + j] - g) * L / (2 * np.pi)) + np.log(np.abs(z[i0 - j] - g) * L / (2 * np.pi))
        far = y - near
        C = np.cov(np.vstack([near, far]))
        vn.append(C[0, 0]); vf.append(C[1, 1]); cv.append(2 * C[0, 1]); vy.append(np.var(y, ddof=1))
    out = {}
    for nm, v in [('Var y', vy), ('Var near', vn), ('Var far', vf), ('2Cov', cv)]:
        v = np.array(v); b = ols(X, v)[nn:]
        rng = np.random.default_rng(5); B = []
        for _ in range(500):
            ix = rng.integers(0, len(rows), len(rows)); Xb = X[ix].copy(); Xb[:, 1:] -= Xb[:, 1:].mean(0)
            B.append(ols(Xb, v[ix])[nn:])
        se = np.std(B, 0, ddof=1)
        print('k=%d %-9s mean %.4f  contrast p=2,3,5,7: ' % (k, nm, np.mean(v)) + ' '.join('%+.4f(%.4f)' % (a, s) for a, s in zip(b[:4], se[:4])), flush=True)
