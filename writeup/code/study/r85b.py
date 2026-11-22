"""Round 85/86 additions after the twenty-fourth review: (1) the Pearson chi^2 of each last-digit table (3 free entries
per row); (2) the uniform 1/4 comparison promised by the registration; (3) the exact first-order chain on classes mod
30030 (one-gap law per class, no memory beyond the class) for the last-digit table at 1e9, k = 1..3, D <= 240.
usage: r85b.py"""
import os, json, math
import numpy as np
import common as C
import r47

OUT = r47.OUT
DIG = (1, 3, 7, 9)
res = {'pearson': {}, 'uniform_max_abs_z': {}}
for name in ('r85_1e9', 'r85_1e12', 'r86_1e15'):
    r = json.load(open(os.path.join(OUT, name + '.json'))); x2 = 0.0; uz = {}
    for k, v in r['k'].items():
        law = np.array(v['law']); dat = np.array(v['data']); N = np.array(v['row_counts'])[:, None]
        x2 += float(np.sum(N * (dat - law) ** 2 / law))
        uz[k] = float(np.max(np.abs((dat - 0.25) / np.sqrt(dat * (1 - dat) / N))))
    res['pearson'][name] = {'chi2': x2, 'dof': 36}; res['uniform_max_abs_z'][name] = uz
# exact first-order class chain at 1e9
PP = r47.PP; CLS = np.array([a for a in range(PP) if math.gcd(a, PP) == 1]); IDX = np.full(PP, -1); IDX[CLS] = np.arange(len(CLS))
LAM = math.log(1.05e9); RHO0 = 1 / (LAM * r47.RHO_DEN); G = 240; NG = G // 2
P1 = np.zeros((len(CLS), NG))
for g in range(2, G + 1, 2):
    st = r47.setup(g); n = st['M'].shape[1]; rows = IDX[CLS[np.gcd(CLS + g, PP) == 1]]; M = st['M']; p = RHO0 * st['s1']
    L = M @ np.log1p(-p) if n else np.zeros(len(rows))
    if n > 1:
        w = p / (1 - p); L = L + 0.5 * np.sum((M @ ((w[:, None] * w[None, :]) * st['R'])) * M, axis=1)
    P1[rows, g // 2 - 1] = RHO0 * st['S2'] * np.exp(L)
P1 /= P1.sum(1, keepdims=True)
r = json.load(open(os.path.join(OUT, 'r85_1e9.json')))
da = np.searchsorted(DIG, CLS % 10); cur = P1.copy(); res['chain_1e9'] = {}
for k in (1, 2, 3):
    if k > 1:
        new = np.zeros_like(cur)
        for j in range(NG):
            nxt = IDX[(CLS + 2 * (j + 1)) % PP]; ok = (nxt >= 0) & (cur[:, j] > 0); m = NG - j - 1
            if m > 0:
                new[ok, j + 1:] += cur[ok, j:j + 1] * P1[nxt[ok], :m]
        cur = new
    T = np.zeros((4, 4))
    for j in range(NG):
        db = np.searchsorted(DIG, (CLS + 2 * (j + 1)) % 10); ok = cur[:, j] > 0
        np.add.at(T, (da[ok], db[ok]), cur[ok, j])
    law = T / T.sum(1, keepdims=True); v = r['k'][str(k)]; dat = np.array(v['data']); N = np.array(v['row_counts'])[:, None]
    z = (law - dat) / np.sqrt(dat * (1 - dat) / N)
    res['chain_1e9'][str(k)] = {'max_abs_z': float(np.max(np.abs(z))), 'sum_z2': float(np.sum(z * z))}
print(json.dumps(res, indent=1))
json.dump(res, open(os.path.join(OUT, 'r85b.json'), 'w'), indent=1)
