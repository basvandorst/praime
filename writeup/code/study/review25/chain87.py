"""Round 87 power check: first-order chain on classes mod 30030 (per-class gap law, round 47/81 z^0 term at
lam = ln(1e18 + 1e6)), repeat probabilities k = 1..6, q = 3..23, D <= 480, against r87.json's law and data."""
import sys, math, json, time
sys.path.insert(0, '/home/bas/praime/writeup/code/study')
import numpy as np, r47
t0 = time.time()
LAM = math.log(1e18 + 1e6); RHO0 = 1 / (LAM * r47.RHO_DEN); DMAX = 480; ND = DMAX // 2
PP = r47.PP; CLS = np.array([a for a in range(PP) if math.gcd(a, PP) == 1]); idx = np.full(PP, -1); idx[CLS] = np.arange(len(CLS))
F = np.zeros((len(CLS), ND)); Ds = np.arange(2, DMAX + 1, 2)
for j, D in enumerate(Ds):
    st = r47.setup(int(D)); M = st['M']; n = M.shape[1]; p = RHO0 * st['s1']
    rows = idx[CLS[np.gcd(CLS + D, PP) == 1]]
    L = M @ np.log1p(-p) if n else np.zeros(len(rows))
    if n > 1:
        w = p / (1 - p); L = L + 0.5 * np.sum((M @ ((w[:, None] * w[None, :]) * st['R'])) * M, axis=1)
    F[rows, j] = RHO0 * st['S2'] * np.exp(L)
print('table', round(time.time() - t0), 's', F.sum(1).min(), F.sum(1).max(), flush=True)
nidx = idx[(CLS[:, None] + Ds[None, :]) % PP]
E = np.zeros((len(CLS), ND))
for j in range(ND):
    m = nidx[:, j] >= 0; np.add.at(E[:, j], nidx[m, j], F[m, j] / len(CLS))
J = json.load(open('/home/bas/praime/writeup/data/study/r87.json')); QS = (3, 5, 7, 11, 13, 17, 19, 23)
tot = 0
for k in range(1, 7):
    dk = E.sum(0); line = []
    for q in QS:
        ch = dk[Ds % q == 0].sum() / dk.sum(); v = J['data'][str(k)][str(q)]
        dl = (J['theory'][str(k)][str(q)] - ch) / v['se']; dd = (v['p'] - ch) / v['se']
        if k >= 2: tot += dd * dd
        line.append('q%d law-chain %+.2f data-chain %+.2f' % (q, dl, dd))
    print(k, '; '.join(line), flush=True)
    if k < 6:
        E2 = np.zeros_like(E)
        for j in range(ND):
            if not E[:, j].any(): continue
            for j2 in range(ND - j - 1):
                w = E[:, j] * F[:, j2]; m = (w > 0) & (nidx[:, j2] >= 0)
                np.add.at(E2[:, j + j2 + 1], nidx[m, j2], w[m])
        E = E2
print('sum (data-chain)^2/se^2 over k>=2 (40 values): %.1f' % tot, round(time.time() - t0), 's')
