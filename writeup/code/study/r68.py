"""Round 68 (preregister/round68.txt): mutual information of consecutive gaps in [1e12, 1e12 + 1e10) against the
discrete leading-order law and I_inf. usage: r68.py"""
import os, json, math
import numpy as np
import common as C
import theory2 as T2

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
L0, L1 = 10 ** 12, 10 ** 12 + 10 ** 10
SH = 4096


def mi_counts(keys, w):
    g1 = keys // SH; g2 = keys % SH
    Hm = lambda c: float(-np.sum((c / c.sum()) * np.log2(c / c.sum())))
    def marg(k):
        u, inv = np.unique(k, return_inverse=True); return np.bincount(inv, weights=w)
    return Hm(marg(g1)) + Hm(marg(g2)) - Hm(w), len(w)


def data():
    acc = {}; prev = None
    for s in range(L0, L1, 10 ** 8):
        P = np.nonzero(C.primes_mask(s, 10 ** 8))[0].astype(np.int64) + s
        if prev is not None:
            P = np.concatenate([prev, P])
        g = np.diff(P); k = g[:-1] * SH + g[1:]
        v, c = np.unique(k, return_counts=True)
        for a, b in zip(v.tolist(), c.tolist()):
            acc[a] = acc.get(a, 0) + b
        prev = P[-2:]
    keys = np.array(list(acc.keys()), np.int64); cnt = np.array(list(acc.values()), float)
    I, cells = mi_counts(keys, cnt); N = cnt.sum()
    nz1 = len(np.unique(keys // SH)); nz2 = len(np.unique(keys % SH))
    I += ((nz1 - 1) + (nz2 - 1) - (cells - 1)) / (2 * N * math.log(2))
    return I, int(N), cells


def law(lam):
    G = int(30 * lam) // 2 * 2
    ev = np.arange(2, G + 1, 2)
    A, B = np.meshgrid(ev, ev, indexing='ij'); A = A.ravel(); B = B.ravel(); m = A + B <= G; A, B = A[m], B[m]
    sets = np.stack([np.zeros_like(A), A, A + B], 1)
    S = np.concatenate([T2.S_many(sets[i:i + 200000]) for i in range(0, len(sets), 200000)])
    w = S * np.exp(-(A + B) / lam); k = w > 0
    I, _ = mi_counts((A * SH + B)[k], w[k])
    return I


if __name__ == '__main__':
    lam = float(np.mean(np.log(np.linspace(L0, L1, 1001))))
    Iinf = json.load(open(os.path.join(OUT, 'r56.json')))['I_inf_bits']
    Il = law(lam); Id, N, cells = data()
    out = {'lambda': lam, 'pairs': N, 'cells': cells, 'I_data': Id, 'I_law_discrete': Il, 'I_inf': Iinf,
           'data_minus_law': Id - Il, 'lambda_times_data_minus_Iinf': lam * (Id - Iinf)}
    print({k: (round(v, 5) if isinstance(v, float) else v) for k, v in out.items()})
    json.dump(out, open(os.path.join(OUT, 'r68.json'), 'w'), indent=1)
