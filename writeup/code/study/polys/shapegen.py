"""Jitter the best optimised root sets (sigma per root), score in batches (same heuristic as shapeopt.py), keep shapes
with score >= SMIN, write 'b4 b3 b2 b1 V' lines (V = 1.3 x max |N|/K on |y| <= 18). usage: shapegen.py IN K N SIGMA SMIN SEED OUT"""
import sys, json, math
import numpy as np
IN, K, NTOT, SIG, SMIN, SEED, OUT = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), int(sys.argv[6]), sys.argv[7]
B5 = 3 if K == 72 else 0
rng = np.random.default_rng(SEED)
base = [np.array([complex(a, b) for a, b in d['roots']]) for d in json.load(open(IN))[:60]]
g = 0.5772156649; kap = math.exp(g) * math.log(67); P2 = 67 * 67; L = 56
yy = np.arange(-70, 71).astype(float); core = np.abs(yy) <= 18
seen = set(); out = open(OUT, 'w'); kept = 0; B = 20000
for start in range(0, NTOT, B):
    idx = rng.integers(0, len(base), B)
    R = np.array([base[i] for i in idx])
    real = np.abs(R.imag) < 1e-9
    jit = rng.normal(0, SIG, R.shape)
    R = R + jit                       # real roots stay real; complex pairs need the same shift for both members
    # restore conjugate pairs: for complex roots, pair k and k+1 (sorted as in shapeopt: pairs at the end)
    for j in range(R.shape[1] - 1):
        pair = (~real[:, j]) & (~real[:, j + 1]) & (np.abs(np.array([base[i][j].real - base[i][j + 1].real for i in idx])) < 1e-9)
        R[pair, j + 1] = np.conj(R[pair, j])
    R = R - ((R.sum(axis=1) + B5) / 6)[:, None]
    C = np.array([np.real(np.poly(r)) for r in R])       # C[:,0] = 1, C[:,1] = B5
    Bi = np.rint(C[:, 2:6]).astype(np.int64)
    V = np.ones((B, len(yy)))
    for j in range(6): V = V * (yy[None, :] - R[:, j][:, None])
    V = np.abs(V.real) / K
    lp = np.where(V < P2, 0.0, np.log(np.minimum(1, kap / np.log(np.maximum(V, 3)))))
    c = np.concatenate([np.zeros((B, 1)), np.cumsum(lp, axis=1)], axis=1)
    sc = (c[:, L:] - c[:, :-L]).max(axis=1) / math.log(10)
    vmax = V[:, core].max(axis=1)
    for i in np.nonzero(sc >= SMIN)[0]:
        t = tuple(Bi[i])
        if t in seen: continue
        seen.add(t); out.write('%d %d %d %d %d\n' % (t + (int(vmax[i] * 1.3) + 1000,))); kept += 1
print('kept', kept, 'of', NTOT)
