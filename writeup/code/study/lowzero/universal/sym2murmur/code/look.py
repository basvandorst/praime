import sys, glob, numpy as np
outdir, k, N0, N1, w = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5])
ys, vs, ws = [], [], []
for f in sorted(glob.glob("%s/data_k%d_N*.npz" % (outdir, k))):
    d = np.load(f); N = int(d["N"])
    if not (N0 <= N <= N1): continue
    p = d["p"].astype(float); dim = int(d["dim"])
    lam2 = d["s"] / dim - 1.0                      # E_nat[lambda_f(p^2)] at level N
    ys.append(p / N ** 2); vs.append(np.sqrt(p) * lam2); ws.append(np.full(len(p), dim, float))
y = np.concatenate(ys); v = np.concatenate(vs); wt = np.concatenate(ws)
edges = np.arange(0, 3.0 + 1e-9, w)
b = np.digitize(y, edges) - 1
print("levels", len(ys), "points", len(y))
for i in range(len(edges) - 1):
    m = b == i
    if m.sum() < 3: continue
    mu = np.sum(wt[m] * v[m]) / np.sum(wt[m])
    # standard error from weighted scatter
    var = np.sum(wt[m] ** 2 * (v[m] - mu) ** 2) / np.sum(wt[m]) ** 2
    print("%5.2f-%5.2f n=%5d  mean %+7.3f  se %.3f" % (edges[i], edges[i + 1], m.sum(), mu, np.sqrt(var)))
