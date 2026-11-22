"""Order 3 of the conditioned expansion against the exact sums (r28.json, g = 74..90)."""
import json, math, os, sys
from itertools import combinations
import numpy as np
import r47
ref = json.load(open(os.path.join(r47.OUT, 'r28.json')))
out = {}
for hs in ref:
    g = int(hs); st = r47.setup(g)
    pts = np.arange(2, g, 2, dtype=np.int64); n = len(pts)
    tri = np.array(list(combinations(range(n), 3)), dtype=np.int64); i, j, k = tri.T; z = np.zeros(len(tri), np.int64)
    s3 = r47.S_gt(np.stack([z, np.full(len(tri), g), pts[i], pts[j], pts[k]], 1)) / st['S2']
    s1, R, M = st['s1'], st['R'], st['M']
    kap = s3 / (s1[i] * s1[j] * s1[k]) - 1 - R[i, j] - R[i, k] - R[j, k]
    for lam in (17.8, 23.0):
        rho0 = 1 / (lam * r47.RHO_DEN); p = rho0 * s1; w = p / (1 - p)
        L1 = M @ np.log1p(-p)
        L2 = L1 + 0.5 * np.sum((M @ (w[:, None] * w[None, :] * R)) * M, axis=1)
        W = M * w[None, :]
        o3 = np.zeros(M.shape[0])
        for c0 in range(0, len(tri), 20000):
            sl = slice(c0, c0 + 20000); o3 += (W[:, i[sl]] * W[:, j[sl]] * W[:, k[sl]]) @ kap[sl]
        L3 = L2 - o3
        C = w[:, None] * w[None, :] * R                      # c~_ij (zero diagonal)
        v = M @ C; v2 = M @ (C * C)
        ov = -0.5 * (0.5 * np.sum(v2 * M, axis=1) + np.sum(M * (v * v - v2), axis=1))
        L2o = L2 + ov; L3o = L3 + ov
        ex = ref[hs]['lam%.1f' % lam]['exact']
        f = lambda L: float(np.sum(rho0 * st['S2'] * np.exp(L)) / st['nclass_all']) / ex - 1
        out['%s_%.1f' % (hs, lam)] = {'o2_rel': f(L2), 'o3_rel': f(L3), 'o2_overlap_rel': f(L2o), 'o3_overlap_rel': f(L3o)}
        print(hs, lam, {kk: round(v, 5) for kk, v in out['%s_%.1f' % (hs, lam)].items()}, flush=True)
json.dump(out, open(os.path.join(r47.OUT, 'r47_check3_y%d.json' % r47.Y), 'w'), indent=1)
