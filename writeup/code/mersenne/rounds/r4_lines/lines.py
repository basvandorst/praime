"""Experiment 2: collinear exponents in the square layout of width w.

n sits at row (n-1)//w, column (n-1)%w.  For every pair i<j of points of a
set the direction (dr, dc) = (r_j - r_i, c_j - c_i) is reduced by
g = gcd(dr, |dc|) (dr >= 0, and dc > 0 when dr = 0, so the reduced pair is a
canonical line direction).  Group the pairs by (set, i, direction): a line
with m points gives its lowest point a group of m-1, its k-th point a group
of m-k, so
   sum over groups of C(size, 2) = sum over lines C(m, 3) = T (collinear triples)
   sum over groups of C(size, 3) = sum over lines C(m, 4) = Q (collinear 4-sets)
   max group size + 1            = Mx (most points on one line)
Per-w p (exchangeable rank over the N+1 sets) and z-scores; family-wise by
min p and by max z over all w and the three statistics.
Usage: lines.py <view> <Nnull> <WMAX>
"""
import sys, os, json, time
import numpy as np
from common import EXPS, HERE

view = sys.argv[1]
NN = int(sys.argv[2])
WMAX = int(sys.argv[3])
null = np.load(os.path.join(HERE, f"null_{view}.npy"))[:NN]
lim = {"full": 10**9, "page": 500000, "h5": 100000}[view]
real = np.array([e for e in EXPS if e < lim], dtype=np.int64)
K = len(real)
X = np.vstack([real[None, :], null])
S, N = X.shape[0], X.shape[0] - 1
I, J = np.triu_indices(K, 1)
SI = (np.arange(S)[:, None] * K + I[None, :]).ravel()   # set*K + i per pair
SH_B, SH_A = 13, 26
assert X.max() // 2 < 2 ** SH_A and 2 * WMAX < 2 ** SH_B
names = ["T", "Q", "Mx"]
NS = 3
minp = np.ones((NS + 1, S))
maxz = np.full((NS + 1, S), -np.inf)
rows = []
t0 = time.time()
for w in range(2, WMAX + 1):
    r = (X - 1) // w
    c = (X - 1) % w
    dr = (r[:, J] - r[:, I]).ravel()
    dc = (c[:, J] - c[:, I]).ravel()
    g = np.gcd(dr, np.abs(dc))
    a = dr // g
    b = dc // g + w
    key = (SI << (SH_A + SH_B)) | (a << SH_B) | b
    key.sort()
    st = np.empty(len(key), dtype=bool)
    st[0] = True
    st[1:] = key[1:] != key[:-1]
    gi = np.cumsum(st) - 1
    cnt = np.bincount(gi)
    gset = (key[st] >> (SH_A + SH_B)) // K
    cf = cnt.astype(np.float64)
    T = np.bincount(gset, cf * (cf - 1) / 2, minlength=S)
    Q = np.bincount(gset, cf * (cf - 1) * (cf - 2) / 6, minlength=S)
    off = np.searchsorted(gset, np.arange(S), side="left")
    Mx = np.maximum.reduceat(cnt, off) + 1
    St = np.vstack([T, Q, Mx.astype(np.float64)])
    rec = {"w": w, "T": int(T[0]), "Q": int(Q[0]), "Mx": int(Mx[0]),
           "null_mean_T": float(T[1:].mean()), "null_mean_Mx": float(Mx[1:].mean())}
    for k in range(NS):
        s = St[k]
        sa = np.sort(s)
        pall = (S - np.searchsorted(sa, s, side="left")) / S
        rec["p_" + names[k]] = float(pall[0])
        minp[k] = np.minimum(minp[k], pall)
        minp[NS] = np.minimum(minp[NS], pall)
        mu, sd = s[1:].mean(), s[1:].std()
        z = (s - mu) / sd if sd > 0 else np.zeros_like(s)
        rec["z_" + names[k]] = float(z[0])
        maxz[k] = np.maximum(maxz[k], z)
        maxz[NS] = np.maximum(maxz[NS], z)
    rows.append(rec)
    if w % 500 == 0:
        print(f"w={w} {time.time()-t0:.0f}s", flush=True)

out = {"view": view, "N": N, "K": K, "WMAX": WMAX, "summary": {}}
for k in range(NS + 1):
    key = names[k] if k < NS else "ALL"
    ks = [names[k]] if k < NS else names
    best = min(((r["p_" + s], r["w"], s) for r in rows for s in ks))
    bz = max(((r["z_" + s], r["w"], s) for r in rows for s in ks))
    out["summary"][key] = {
        "min_p": best[0], "w": best[1], "stat": best[2],
        "familywise_p_minp": int((minp[k] <= minp[k][0] + 1e-12).sum()) / S,
        "max_z": bz[0], "w_z": bz[1], "stat_z": bz[2],
        "familywise_p_maxz": int((maxz[k] >= maxz[k][0] - 1e-12).sum()) / S}
out["top"] = sorted(rows, key=lambda r: min(r["p_" + s] for s in names))[:15]
out["w708"] = rows[708 - 2] if WMAX >= 708 else None
with open(os.path.join(HERE, f"lines_{view}.json"), "w") as f:
    json.dump(out, f, indent=1)
with open(os.path.join(HERE, f"lines_{view}_perw.json"), "w") as f:
    json.dump(rows, f)
for k, v in out["summary"].items():
    print(k, v)
for r in out["top"]:
    print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})
print("w=708:", out["w708"])
