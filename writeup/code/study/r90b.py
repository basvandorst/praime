"""Round 90 additions after the twenty-sixth review: (1) the law at a single lam equal to the primes' mean ln p in
[1e7, 1e8); (2) the exact first-order class chain mod 30030 at that lam, as comparator; (3) the q <= 13 / q >= 17 split of
r82c. usage: r90b.py"""
import os, json, math
import numpy as np
import common as C
import r47

OUT = r47.OUT; PP = r47.PP
CLS = np.array([a for a in range(PP) if math.gcd(a, PP) == 1]); IDX = np.full(PP, -1); IDX[CLS] = np.arange(len(CLS))
DMAX = 200; NZ = 48; KS = (1, 2, 3, 4); QS = (3, 5, 7, 11, 13)
Z = np.exp(2j * np.pi * np.arange(NZ) / NZ)
P = np.nonzero(C.primes_mask(10 ** 7, 9 * 10 ** 7))[0].astype(np.int64) + 10 ** 7
LAM = float(np.mean(np.log(P.astype(float)))); RHO0 = 1 / (LAM * r47.RHO_DEN)
r90 = json.load(open(os.path.join(OUT, 'r90.json')))
law = {k: np.zeros(DMAX // 2) for k in KS}; P1 = np.zeros((len(CLS), DMAX // 2))
for D in range(2, DMAX + 1, 2):
    st = r47.setup(D); M = st['M']; n = M.shape[1]; p = RHO0 * st['s1']
    logG = np.zeros((M.shape[0], NZ), complex)
    if n:
        logG += M @ np.log1p(p[:, None] * (Z[None, :] - 1))
    if n > 1:
        W = p[:, None] * (Z[None, :] - 1) / (1 + p[:, None] * (Z[None, :] - 1))
        for m in range(NZ):
            X = M * W[:, m][None, :]
            logG[:, m] += 0.5 * np.sum((X.real @ st['R'] + 1j * (X.imag @ st['R'])) * X, axis=1)
    coef = np.fft.fft(np.exp(logG), axis=1).real / NZ
    for k in KS:
        law[k][D // 2 - 1] = RHO0 * st['S2'] * coef[:, k - 1].sum()
    P1[IDX[CLS[np.gcd(CLS + D, PP) == 1]], D // 2 - 1] = RHO0 * st['S2'] * coef[:, 0]
Ds = np.arange(2, DMAX + 1, 2); res = {'mean_ln_p': LAM, 'law_single_lam': {}, 'chain': {}}
cur = P1.copy(); zl = []; zc = []
for k in KS:
    if k > 1:
        new = np.zeros_like(cur)
        for j in range(DMAX // 2):
            nxt = IDX[(CLS + 2 * (j + 1)) % PP]; ok = (nxt >= 0) & (cur[:, j] > 0); m = DMAX // 2 - j - 1
            if m > 0:
                new[ok, j + 1:] += cur[ok, j:j + 1] * P1[nxt[ok], :m]
        cur = new
    ch = cur.mean(0); res['law_single_lam'][str(k)] = {}; res['chain'][str(k)] = {}
    for q in QS:
        v = r90['data'][str(k)][str(q)]
        a = float(law[k][Ds % q == 0].sum() / law[k].sum()); c = float(ch[Ds % q == 0].sum() / ch.sum())
        res['law_single_lam'][str(k)][str(q)] = {'p': a, 'z': (a - v['p']) / v['se']}
        res['chain'][str(k)][str(q)] = {'p': c, 'z': (c - v['p']) / v['se']}
        zl.append((a - v['p']) / v['se'])
        if k >= 2:
            zc.append((c - v['p']) / v['se'])
res['law_single_lam_max_abs_z'] = float(np.max(np.abs(zl))); res['law_single_lam_sum_z2'] = float(np.sum(np.square(zl)))
res['chain_k_ge_2_max_abs_z'] = float(np.max(np.abs(zc))); res['chain_k_ge_2_sum_z2'] = float(np.sum(np.square(zc)))
# split of r82c
r = json.load(open(os.path.join(OUT, 'r82c.json'))); res['r82c_split'] = {}
for name in ('r82', 'r84'):
    for lab, qs in (('q<=13', ('3', '5', '7', '11', '13')), ('q>=17', ('17', '19', '23'))):
        zlaw = [r[name]['per_k'][k][q]['z_law'] for k in map(str, range(2, 7)) for q in qs]
        zch = [r[name]['per_k'][k][q]['z_chain'] for k in map(str, range(2, 7)) for q in qs]
        res['r82c_split'][name + ' ' + lab] = {'n': len(zlaw), 'law_sum_z2': float(np.sum(np.square(zlaw))), 'law_max': float(np.max(np.abs(zlaw))),
                                               'chain_sum_z2': float(np.sum(np.square(zch))), 'chain_max': float(np.max(np.abs(zch)))}
print({k: v for k, v in res.items() if k not in ('law_single_lam', 'chain')})
json.dump(res, open(os.path.join(OUT, 'r90b.json'), 'w'), indent=1)
