"""Round 48 (preregister/round48.txt): the counts of large gaps (120..298) in [1e9, 2e9) against the conditioned
expansion of round 47 (order 3 + overlap correction, y = 13). usage: r48.py"""
import os, json, math
import numpy as np
import mpmath
import common as C
import r47

L0, L1 = 10 ** 9, 2 * 10 ** 9
LAMS = np.array([math.log(L0), math.log(1.5e9), math.log(L1)])
GS = list(range(120, 299, 2))
BANDS = [(120, 178), (180, 238), (240, 298)]


def shares3(g):
    """ln share at LAMS with order 2, and with order 3 + overlap correction."""
    st = r47.setup(g)
    pts = np.arange(2, g, 2, dtype=np.int64); n = len(pts)
    s1, R, M = st['s1'], st['R'], st['M']
    K = np.zeros((n, n, n))
    for i in range(n - 2):
        jj, kk = np.triu_indices(n - i - 1, 1); jj = jj + i + 1; kk = kk + i + 1
        z = np.zeros(len(jj), np.int64)
        s3 = r47.S_gt(np.stack([z, np.full(len(jj), g), np.full(len(jj), pts[i]), pts[jj], pts[kk]], 1)) / st['S2']
        kap = s3 / (s1[i] * s1[jj] * s1[kk]) - 1 - R[i, jj] - R[i, kk] - R[jj, kk]
        K[i, jj, kk] = kap; K[i, kk, jj] = kap
    o2, o3 = [], []
    for lam in LAMS:
        rho0 = 1 / (lam * r47.RHO_DEN); p = rho0 * s1; w = p / (1 - p)
        Lb = M @ np.log1p(-p)
        Cm = w[:, None] * w[None, :] * R
        L2 = Lb + 0.5 * np.sum((M @ Cm) * M, axis=1)
        W = M * w[None, :]
        t3 = np.zeros(M.shape[0])
        for i in range(n - 2):
            Wi = W[:, i + 1:]
            t3 += W[:, i] * 0.5 * np.sum((Wi @ K[i, i + 1:, i + 1:]) * Wi, axis=1)
        v = M @ Cm; v2 = M @ (Cm * Cm)
        ov = -0.5 * (0.5 * np.sum(v2 * M, axis=1) + np.sum(M * (v * v - v2), axis=1))
        f = rho0 * st['S2'] / st['nclass_all']
        o2.append(math.log(f * np.sum(np.exp(L2)))); o3.append(math.log(f * np.sum(np.exp(L2 - t3 + ov))))
    return o2, o3


def counts():
    cnt = {}; last = None; ngaps = 0
    for s in range(L0, L1, 10 ** 7):
        P = np.nonzero(C.primes_mask(s, 10 ** 7))[0].astype(np.int64) + s
        if last is not None:
            P = np.concatenate([[last], P])
        last = int(P[-1])
        g = np.diff(P); ngaps += len(g)
        v, c = np.unique(g[g >= 120], return_counts=True)
        for a, b in zip(v.tolist(), c.tolist()):
            cnt[a] = cnt.get(a, 0) + b
    return cnt, ngaps


if __name__ == '__main__':
    cnt, ngaps = counts()
    norm = ngaps / float(mpmath.li(L1) - mpmath.li(L0))
    u = np.linspace(math.log(L0), math.log(L1), 2001)
    per = {}
    for g in GS:
        o2, o3 = shares3(g)
        e = {}
        for name, ls in (('o2', o2), ('o3', o3)):
            sh = np.exp(np.interp(u, LAMS, ls))
            e[name] = float(np.trapezoid(sh / u * np.exp(u), u)) * norm
        per[str(g)] = {'count': cnt.get(g, 0), 'law_o2': e['o2'], 'law_o3': e['o3']}
        print(g, per[str(g)], flush=True)
    res = {'gaps': ngaps, 'norm': norm, 'per_gap': per, 'bands': {}, 'by_class': {}}
    for lo, hi in BANDS:
        c = sum(per[str(g)]['count'] for g in GS if lo <= g <= hi)
        for name in ('o2', 'o3'):
            lw = sum(per[str(g)]['law_' + name] for g in GS if lo <= g <= hi)
            res['bands']['%d-%d_%s' % (lo, hi, name)] = {'count': c, 'law': lw, 'ratio': c / lw, 'z': (c - lw) / math.sqrt(lw)}
    for lo, hi in BANDS:
        for cls, f in (('six', lambda g: g % 6 == 0), ('not_six', lambda g: g % 6 != 0)):
            c = sum(per[str(g)]['count'] for g in GS if lo <= g <= hi and f(g))
            lw = sum(per[str(g)]['law_o3'] for g in GS if lo <= g <= hi and f(g))
            res['by_class']['%d-%d_%s' % (lo, hi, cls)] = {'count': c, 'law': lw, 'ratio': c / lw, 'z': (c - lw) / math.sqrt(lw)}
    print(json.dumps(res['bands'], indent=0), json.dumps(res['by_class'], indent=0), flush=True)
    json.dump(res, open(os.path.join(r47.OUT, 'r48.json'), 'w'), indent=1)
