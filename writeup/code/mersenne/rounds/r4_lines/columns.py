"""Experiment 1: column concentration (residues mod w), w = 2..WMAX.

For each set (row 0 = real exponents, rows 1..N = null) and each w:
  D   = number of distinct columns            (small = concentrated)
  M   = largest column count                  (large = concentrated)
  CHI = chi^2 against uniform over the phi(w) coprime columns, on the
        n_c elements coprime to w (a prime dividing w sits alone in its
        column and is left out of chi^2 and G)
  G   = likelihood-ratio statistic 2 sum n_r ln(n_r / E), same classes
Per-w empirical p: real p = (1 + #{null >= real}) / (N + 1); null i p =
#{null j >= null i} / N.  Family-wise p = share of null sets whose min p
(over all w and the chosen statistics) is <= the real min p.
Usage: columns.py <view> [WMAX]   view in full, page, h5
"""
import sys, os, json
import numpy as np
from sympy import totient
from common import EXPS, HERE

view = sys.argv[1]
WMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 10000
null = np.load(os.path.join(HERE, f"null_{view}.npy"))
lim = {"full": 10**9, "page": 500000, "h5": 100000}[view]
real = np.array([e for e in EXPS if e < lim], dtype=np.int64)
K = len(real)
assert null.shape[1] == K
X = np.vstack([real[None, :], null])          # (N+1, K)
N = null.shape[0]
NS = 4
names = ["D_low", "M_high", "CHI_high", "G_high"]

minp = np.ones((NS + 1, N + 1))     # per stat, and overall (last); col 0 = real
maxz = np.full((NS + 1, N + 1), -np.inf)
rows = []                                      # per w: real stats and p's
for w in range(2, WMAX + 1):
    R = np.sort(X % w, axis=1)
    st = np.ones_like(R, dtype=bool)
    st[:, 1:] = R[:, 1:] != R[:, :-1]
    runid = np.cumsum(st.ravel()) - 1
    cnt = np.bincount(runid)
    D = st.sum(axis=1)
    off = np.concatenate([[0], np.cumsum(D)[:-1]])
    M = np.maximum.reduceat(cnt, off)
    sumsq = np.add.reduceat(cnt * cnt, off)
    cf = cnt.astype(np.float64)
    nln = np.add.reduceat(cf * np.log(cf), off)
    nonc = (w % X == 0).sum(axis=1)            # primes dividing w
    nc = K - nonc
    phi = int(totient(w))
    CHI = phi / nc * (sumsq - nonc) - nc
    G = 2.0 * (nln - nc * np.log(nc / phi))
    S = np.vstack([-D, M, CHI, G]).astype(np.float64)   # larger = more striking
    S = np.round(S, 9)
    rec = {"w": w, "D": int(D[0]), "M": int(M[0]), "CHI": float(CHI[0]),
           "G": float(G[0])}
    for k in range(NS):
        s = S[k]
        sa = np.sort(s)                       # all N+1 sets, exchangeable
        ge = (N + 1) - np.searchsorted(sa, s, side="left")
        pall = ge / (N + 1)
        rec["p_" + names[k]] = float(pall[0])
        minp[k] = np.minimum(minp[k], pall)
        minp[NS] = np.minimum(minp[NS], pall)
        mu, sd = s[1:].mean(), s[1:].std()
        z = (s - mu) / sd if sd > 0 else np.zeros_like(s)
        rec["z_" + names[k]] = float(z[0])
        maxz[k] = np.maximum(maxz[k], z)
        maxz[NS] = np.maximum(maxz[NS], z)
    rows.append(rec)

out = {"view": view, "N": N, "K": K, "WMAX": WMAX, "summary": {}}
for k in range(NS + 1):
    key = names[k] if k < NS else "ALL"
    ks = [names[k]] if k < NS else names
    best = min(((r["p_" + s], r["w"], s) for r in rows for s in ks))
    bz = max(((r["z_" + s], r["w"], s) for r in rows for s in ks))
    fw = int((minp[k] <= minp[k][0] + 1e-12).sum()) / (N + 1)
    fwz = int((maxz[k] >= maxz[k][0] - 1e-12).sum()) / (N + 1)
    out["summary"][key] = {"min_p": best[0], "w": best[1], "stat": best[2],
                           "familywise_p_minp": fw,
                           "share_null_reaching_min_p": float((minp[k][1:] <= minp[k][0] + 1e-12).mean()),
                           "max_z": bz[0], "w_z": bz[1], "stat_z": bz[2],
                           "familywise_p_maxz": fwz}
# top 15 w by overall min p
top = sorted(rows, key=lambda r: -max(r["z_" + s] for s in names))[:15]
out["top"] = top
out["w708"] = rows[708 - 2] if WMAX >= 708 else None
with open(os.path.join(HERE, f"columns_{view}.json"), "w") as f:
    json.dump(out, f, indent=1)
np.savez(os.path.join(HERE, f"columns_{view}_fw.npz"), minp=minp, maxz=maxz)
with open(os.path.join(HERE, f"columns_{view}_perw.json"), "w") as f:
    json.dump(rows, f)
print(json.dumps(out["summary"], indent=1))
for r in top:
    print({k: (round(v, 5) if isinstance(v, float) else v) for k, v in r.items()})
print("w=708:", out["w708"])
