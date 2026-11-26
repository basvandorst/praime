"""(1) empirical bin errors vs the independent-forms prediction sigma^2 = sum dim p / (sum dim)^2;
   (2) correlation of k=2 and k=4 residuals at the same (N, p)."""
import sys, glob
sys.path.insert(0, sys.argv[1])
import numpy as np
from model import M_nat
outdir = sys.argv[2]
grid = np.linspace(1e-5, 3.01, 30001)
for k, N0, N1 in ((2, 100, 200), (4, 100, 170)):
    g = M_nat(grid, k)
    ys, rs, ws, ps = [], [], [], []
    for f in sorted(glob.glob("%s/data_k%d_N*.npz" % (outdir, k))):
        d = np.load(f); N = int(d["N"])
        if not (N0 <= N <= N1): continue
        p = d["p"].astype(float); dim = int(d["dim"]); y = p / N ** 2
        rho = (k - 1) * N / (12.0 * dim)
        r = np.sqrt(p) * (d["s"] / dim - 1.0) - rho * (1 / np.sqrt(p) + np.interp(y, grid, g))
        ys.append(y); rs.append(r); ws.append(np.full(len(p), float(dim))); ps.append(p)
    y = np.concatenate(ys); r = np.concatenate(rs); w = np.concatenate(ws); p = np.concatenate(ps)
    # normalized residual z = r / sqrt(p/dim)
    z = r / np.sqrt(p / w)
    print("k=%d: residual / sqrt(p/dim): mean %.4f, var %.4f (independent Sato-Tate forms: 1)" % (k, z.mean(), z.var()))
    for lo, hi in ((0, 0.5), (0.5, 1.5), (1.5, 3.0)):
        m = (y > lo) & (y <= hi)
        print("    y in (%.1f,%.1f]: var %.4f" % (lo, hi, z[m].var()))
# correlation k=2 vs k=4
c = []
for N in range(100, 171):
    f2, f4 = "%s/data_k2_N%d.npz" % (outdir, N), "%s/data_k4_N%d.npz" % (outdir, N)
    try:
        a, b = np.load(f2), np.load(f4)
    except FileNotFoundError:
        continue
    za = a["s"] / int(a["dim"]) - 1.0; zb = b["s"] / int(b["dim"]) - 1.0
    c.append(np.corrcoef(za, zb)[0, 1])
print("corr(k=2, k=4) of E[lambda(p^2)] at the same (N,p), per level: mean %.3f, range %.3f..%.3f" % (np.mean(c), min(c), max(c)))
