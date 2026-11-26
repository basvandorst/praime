"""Exact decomposition of the data (review). For p != N: sum_f lambda_f(p^2) = Tr T_{p^2}/p^{k-1}, and (Cohen/Eichler-Selberg,
level N minus level 1) Tr T_{p^2} = p^{k-2}(k-1)N/12 - (2+p^{k-1})/2 + A2chi + A2N2,
A2N2 = -N sum_{1<=m<4p/N^2} p^{k-2} U_{k-2}(1-m/(2y)) H(m(4p-mN^2)).
Normalized V = sqrt(p)(s/dim-1) = plancherel + hyperbolic + CN_exact + R, R = normalized A2chi (remainder).
Compares the binned CN_exact with its nu-average and R with 24 rho y^{3/2} [k=2] (0 for k=4); measures where the noise lives."""
import sys, os, json
import numpy as np, cypari2
from scipy.special import eval_chebyu
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
DD = "/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/sym2murmur/out_holdout"
pari = cypari2.Pari(); pari.allocatemem(300 * 10 ** 6)
k = int(sys.argv[1]); levels = [int(x) for x in sys.argv[2].split(",")]
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "nu_only.py")).read())
edges = np.arange(0, 3.0 + 1e-9, 0.05); NB = 60
acc = {q: np.zeros(NB) for q in ("w", "V", "sm", "cnx", "cnp", "R", "Rp", "V2", "R2", "dcn2", "cnx2")}
for N in levels:
    d = np.load(DD + "/data_k%d_N%d.npz" % (k, N)); dim = int(d["dim"]); P = d["p"].astype(np.int64)
    p = P.astype(float); y = p / N ** 2; rho = (k - 1) * N / (12.0 * dim)
    V = np.sqrt(p) * (d["s"] / dim - 1.0)
    sm = rho / np.sqrt(p) - 6 * rho * np.sqrt(y) * (1 + 2 * p ** (1 - k)) / (k - 1)
    cnx = np.zeros_like(p); cnp = np.zeros_like(p)
    for m in range(1, 12):
        msk = (4 * P - m * N * N) > 0
        if not msk.any(): continue
        Dv = m * (4 * P[msk] - m * N * N)
        H = np.array([float(x) for x in pari("[" + ",".join("qfbhclassno(%d)" % x for x in Dv) + "]")])
        U = eval_chebyu(k - 2, 1 - m / (2 * y[msk]))
        cnx[msk] += -(12 * rho / (k - 1)) * U * H / np.sqrt(p[msk])
        cnp[msk] += -(12 * rho / (np.pi * (k - 1) * np.sqrt(y[msk]))) * nu(m) * np.sqrt(m * (4 * y[msk] - m)) * U
    R = V - sm - cnx
    Rp = (24 * rho * y ** 1.5) if k == 2 else 0 * y
    b = np.digitize(y, edges) - 1; ok = (b >= 0) & (b < NB); b = b[ok]; w = np.full(ok.sum(), float(dim))
    for q, arr in (("V", V), ("sm", sm), ("cnx", cnx), ("cnp", cnp), ("R", R), ("Rp", Rp)):
        acc[q] += np.bincount(b, w * arr[ok], NB)
    acc["w"] += np.bincount(b, w, NB)
    acc["V2"] += np.bincount(b, w * w * V[ok] ** 2, NB); acc["R2"] += np.bincount(b, w * w * R[ok] ** 2, NB)
    acc["dcn2"] += np.bincount(b, w * w * (cnx - cnp)[ok] ** 2, NB); acc["cnx2"] += np.bincount(b, w * w * cnx[ok] ** 2, NB)
    print("N=%d done" % N, flush=True)
np.savez(os.path.join(os.path.dirname(os.path.abspath(__file__)), "r159_exact_k%d.npz" % k), levels=levels, **acc)
