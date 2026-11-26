import sys, glob, numpy as np
outdir, k, N0, N1, w, ymax = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6])
def load(cls):
    ys, vs, ws = [], [], []
    for f in sorted(glob.glob("%s/data_k%d_N*.npz" % (outdir, k))):
        d = np.load(f); N = int(d["N"])
        if not (N0 <= N <= N1) or (cls is not None and N % 4 != cls): continue
        p = d["p"].astype(float); dim = int(d["dim"])
        ys.append(p / N ** 2); vs.append(np.sqrt(p) * (d["s"] / dim - 1.0)); ws.append(np.full(len(p), dim, float))
    return np.concatenate(ys), np.concatenate(vs), np.concatenate(ws)
edges = np.arange(0, ymax + 1e-9, w)
res = {}
for cls in (None, 1, 3):
    y, v, wt = load(cls)
    b = np.digitize(y, edges) - 1
    out = []
    for i in range(len(edges) - 1):
        m = b == i
        mu = np.sum(wt[m] * v[m]) / np.sum(wt[m]); se = np.sqrt(np.sum(wt[m] ** 2 * (v[m] - mu) ** 2)) / np.sum(wt[m])
        out.append((mu, se))
    res[cls] = out
print("k=%d N in [%d,%d]   all | N=1 mod 4 | N=3 mod 4" % (k, N0, N1))
for i in range(len(edges) - 1):
    print("%5.3f-%5.3f  %+6.2f(%.2f)  %+6.2f(%.2f)  %+6.2f(%.2f)" % (edges[i], edges[i+1], *res[None][i], *res[1][i], *res[3][i]))
