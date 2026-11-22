"""Comparator added after the twenty-fifth review (lesson: run the strongest fair comparator): the exact first-order
chain on the classes mod 30030 (one-gap law per class at the window's lam, no memory beyond the class) for the repeat
probabilities of rounds 82 (1e12, D <= 320) and 84 (1e15, D <= 400), k = 1..6, q = 3..23, conditioned on
p_{n+k} - p_n <= DMAX (paths within DMAX), against the data and the k-th neighbour law stored in r82.json, r84.json.
usage: r82c.py"""
import os, json, math
import numpy as np
import r47

assert r47.Y == 13
OUT = r47.OUT
PP = r47.PP; CLS = np.array([a for a in range(PP) if math.gcd(a, PP) == 1]); IDX = np.full(PP, -1); IDX[CLS] = np.arange(len(CLS))
QS = (3, 5, 7, 11, 13, 17, 19, 23); KS = tuple(range(1, 7))
res = {}
for name, lam, G in (('r82', math.log(1.001e12), 320), ('r84', math.log(1.0000005e15), 400)):
    RHO0 = 1 / (lam * r47.RHO_DEN); NG = G // 2; P1 = np.zeros((len(CLS), NG))
    for g in range(2, G + 1, 2):
        st = r47.setup(g); n = st['M'].shape[1]; rows = IDX[CLS[np.gcd(CLS + g, PP) == 1]]; M = st['M']; p = RHO0 * st['s1']
        L = M @ np.log1p(-p) if n else np.zeros(len(rows))
        if n > 1:
            w = p / (1 - p); L = L + 0.5 * np.sum((M @ ((w[:, None] * w[None, :]) * st['R'])) * M, axis=1)
        P1[rows, g // 2 - 1] = RHO0 * st['S2'] * np.exp(L)
    Ds = np.arange(2, G + 1, 2); r = json.load(open(os.path.join(OUT, name + '.json'))); cur = P1.copy(); out = {}
    zl, zc = [], []
    for k in KS:
        if k > 1:
            new = np.zeros_like(cur)
            for j in range(NG):
                nxt = IDX[(CLS + 2 * (j + 1)) % PP]; ok = (nxt >= 0) & (cur[:, j] > 0); m = NG - j - 1
                if m > 0:
                    new[ok, j + 1:] += cur[ok, j:j + 1] * P1[nxt[ok], :m]
            cur = new
        dk = cur.mean(0); out[str(k)] = {}
        for q in QS:
            c = float(dk[Ds % q == 0].sum() / dk.sum()); v = r['data'][str(k)][str(q)]
            out[str(k)][str(q)] = {'chain': c, 'z_chain': (c - v['p']) / v['se'], 'z_law': v['z']}
            if k >= 2:
                zl.append(v['z']); zc.append((c - v['p']) / v['se'])
    res[name] = {'lambda': lam, 'DMAX': G, 'per_k': out, 'k_ge_2': {'sum_z2_law': float(np.sum(np.square(zl))), 'sum_z2_chain': float(np.sum(np.square(zc))),
                 'max_abs_z_law': float(np.max(np.abs(zl))), 'max_abs_z_chain': float(np.max(np.abs(zc))), 'values': len(zl)}}
    print(name, res[name]['k_ge_2'], flush=True)
json.dump(res, open(os.path.join(OUT, 'r82c.json'), 'w'), indent=1)
