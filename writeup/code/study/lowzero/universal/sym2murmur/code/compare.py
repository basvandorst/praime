"""Compare natural-weight data E[sqrt(p) lambda_f(p^2)] with model densities, binned in y = p/N^2.
Model at (N,p): rho_N (1/sqrt p + M(y)), rho_N = (k-1) N / (12 dim).  Usage: compare.py CODEDIR OUTDIR K N0 N1 BINW YMAX"""
import sys, glob
sys.path.insert(0, sys.argv[1])
import numpy as np
from model import M_nat, M_harm, K21
from kernels import K3, K2
outdir, k, N0, N1, bw, ymax = sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7])
grid = np.linspace(1e-5, ymax + 0.01, 40001)
models = {
    "M_nat (formula)": M_nat(grid, k),
    "no hyperbolic term (M^P)": M_nat(grid, k, hyperbolic=False),
    "smooth terms only": (24 * grid ** 1.5 if k == 2 else 0 * grid) - 6 * np.sqrt(grid) / (k - 1),
    "null (M = 0)": 0 * grid,
}
g3 = np.linspace(1e-4, ymax + 0.01, 1201)
k3 = np.array([K3(z, k) for z in g3])
models["0.6364*K3 (GL(3) main line only)"] = 0.636446 * np.interp(grid, g3, k3)
models["0.6364*K2 (degree-2 kernel)"] = 0.636446 * K2(grid, k)
ys, ds, ws, rhos, ps = [], [], [], [], []
for f in sorted(glob.glob("%s/data_k%d_N*.npz" % (outdir, k))):
    d = np.load(f); N = int(d["N"])
    if not (N0 <= N <= N1): continue
    p = d["p"].astype(float); dim = int(d["dim"])
    ys.append(p / N ** 2); ds.append(np.sqrt(p) * (d["s"] / dim - 1.0)); ws.append(np.full(len(p), float(dim)))
    rhos.append(np.full(len(p), (k - 1) * N / (12.0 * dim))); ps.append(p)
y = np.concatenate(ys); dat = np.concatenate(ds); wt = np.concatenate(ws); rho = np.concatenate(rhos); p = np.concatenate(ps)
sel = y <= ymax
y, dat, wt, rho, p = y[sel], dat[sel], wt[sel], rho[sel], p[sel]
edges = np.arange(0, ymax + 1e-9, bw); b = np.digitize(y, edges) - 1; nb = len(edges) - 1
def binned(v):
    return np.array([np.sum(wt[b == i] * v[b == i]) / np.sum(wt[b == i]) for i in range(nb)])
D = binned(dat)
se = np.array([np.sqrt(np.sum(wt[b == i] ** 2 * (dat[b == i] - D[i]) ** 2)) / np.sum(wt[b == i]) for i in range(nb)])
print("k=%d, %d levels N in [%d,%d], %d (N,p) pairs, %d bins of width %.3f on (0,%.1f]" % (k, len(ys), N0, N1, len(y), nb, bw, ymax))
for name, g in models.items():
    mv = rho * (1 / np.sqrt(p) + np.interp(y, grid, g))
    Mb = binned(mv)
    chi2 = np.sum(((D - Mb) / se) ** 2)
    a = np.sum(D * Mb / se ** 2) / np.sum(Mb ** 2 / se ** 2) if np.any(Mb != 0) else np.nan
    ae = 1 / np.sqrt(np.sum(Mb ** 2 / se ** 2)) if np.any(Mb != 0) else np.nan
    print("  %-34s chi2 = %7.1f   fitted amplitude %.3f +- %.3f" % (name, chi2, a, ae))
