"""Discovery: two-component fit M = a*smooth + b*classnumber-part; also chi2 for nu=1 and Zubrilina-shape comparators."""
import sys, glob
sys.path.insert(0, sys.argv[1])
import numpy as np
from model import M_parts, M_nat_nu1, M_zub, M_nat
outdir = sys.argv[2]
grid = np.linspace(1e-5, 3.01, 30001)
for k, N0, N1 in ((2, 100, 200), (4, 100, 170)):
    sm, cn = M_parts(grid, k); nu1 = M_nat_nu1(grid, k); zb = M_zub(grid, k)
    ys, ds, ws, rh, ps = [], [], [], [], []
    for f in sorted(glob.glob("%s/data_k%d_N*.npz" % (outdir, k))):
        d = np.load(f); N = int(d["N"])
        if not (N0 <= N <= N1): continue
        p = d["p"].astype(float); dim = int(d["dim"])
        ys.append(p / N ** 2); ds.append(np.sqrt(p) * (d["s"] / dim - 1.0)); ws.append(np.full(len(p), float(dim)))
        rh.append(np.full(len(p), (k - 1) * N / (12.0 * dim))); ps.append(p)
    y, D, w, rho, p = map(np.concatenate, (ys, ds, ws, rh, ps))
    edges = np.arange(0, 3.0 + 1e-9, 0.05); b = np.digitize(y, edges) - 1; nb = len(edges) - 1
    bm = lambda v: np.array([np.sum(w[b == i] * v[b == i]) / np.sum(w[b == i]) for i in range(nb)])
    Db = bm(D); se = np.array([np.sqrt(np.sum(w[b == i] ** 2 * (D[b == i] - Db[i]) ** 2)) / np.sum(w[b == i]) for i in range(nb)])
    pred_se = np.array([np.sqrt(0.96 * np.sum(w[b == i] * p[b == i])) / np.sum(w[b == i]) for i in range(nb)])
    base = bm(rho / np.sqrt(p)); S = bm(rho * np.interp(y, grid, sm)); C = bm(rho * np.interp(y, grid, cn))
    X = np.vstack([S, C]).T / se[:, None]; t = (Db - base) / se
    cov = np.linalg.inv(X.T @ X); ab = cov @ X.T @ t
    print("k=%d: a(smooth) = %.3f +- %.3f, b(class numbers) = %.3f +- %.3f, corr %.2f;  empirical/predicted se: median %.3f" % (
        k, ab[0], np.sqrt(cov[0, 0]), ab[1], np.sqrt(cov[1, 1]), cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1]), np.median(se / pred_se)))
    for name, g in (("formula", sm + cn), ("nu = 1", nu1), ("Zubrilina shape", zb)):
        Mb = base + bm(rho * np.interp(y, grid, g))
        print("     %-16s chi2 = %.1f" % (name, np.sum(((Db - Mb) / se) ** 2)))
