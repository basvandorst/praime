"""Round 58 (preregister/round58.txt): conditional mutual information I(g_{n+1}; g_{n-1} | g_n) in long windows, against
J_inf (r58_theory.json) and against the discrete leading-order law P(a,b,c) ~ S4 e^(-(a+b+c)/lam) at the window's lam.
usage: r58.py data | law"""
import os, sys, json, math
import numpy as np
import common as C
import theory2 as T2
import theory_r3 as R3

OUT = R3.OUT; PATH = os.path.join(OUT, 'r58.json')
WINS = {'train': (10 ** 7, 9 * 10 ** 7), '1e9-2e9': (10 ** 9, 10 ** 9), '1e10-2e10': (10 ** 10, 10 ** 10)}
SH = 2048


def Hmm(c):
    c = np.asarray(c, float); n = c.sum(); p = c / n
    return float(-np.sum(p * np.log2(p)) + (len(c) - 1) / (2 * n * math.log(2)))


def cmi_from_triples(keys, cnt):
    g1 = keys // (SH * SH); g2 = (keys // SH) % SH; g3 = keys % SH
    def marg(k):
        u, inv = np.unique(k, return_inverse=True); return np.bincount(inv, weights=cnt)
    return Hmm(marg(g1 * SH + g2)) + Hmm(marg(g2 * SH + g3)) - Hmm(marg(g2)) - Hmm(cnt)


def data(L, W):
    acc = {}; prev = None; seg = 10 ** 8
    for s in range(L, L + W, seg):
        P = np.nonzero(C.primes_mask(s, min(seg, L + W - s)))[0].astype(np.int64) + s
        if prev is not None:
            P = np.concatenate([prev, P])
        g = np.diff(P)
        k = (g[:-2] * SH + g[1:-1]) * SH + g[2:]
        v, c = np.unique(k, return_counts=True)
        for a, b in zip(v.tolist(), c.tolist()):
            acc[a] = acc.get(a, 0) + b
        prev = P[-3:]
    keys = np.array(list(acc.keys()), dtype=np.int64); cnt = np.array(list(acc.values()), float)
    return cmi_from_triples(keys, cnt), int(cnt.sum()), len(cnt)


def law(lam, G):
    """discrete leading-order law over even a, b, c with a + b + c <= G: P ~ S({0,a,a+b,a+b+c}) e^(-(a+b+c)/lam)."""
    ev = np.arange(2, G + 1, 2)
    A, B, Cc = np.meshgrid(ev, ev, ev, indexing='ij'); A = A.ravel(); B = B.ravel(); Cc = Cc.ravel()
    m = A + B + Cc <= G; A, B, Cc = A[m], B[m], Cc[m]
    sets = np.stack([np.zeros_like(A), A, A + B, A + B + Cc], 1)
    S = np.concatenate([T2.S_many(sets[i:i + 200000]) for i in range(0, len(sets), 200000)])
    w = S * np.exp(-(A + B + Cc) / lam)
    keys = (A * SH + B) * SH + Cc
    k = w > 0
    return cmi_from_triples(keys[k], w[k] * 1e12)            # weights as large pseudo-counts: Miller-Madow term negligible


if __name__ == '__main__':
    res = json.load(open(PATH)) if os.path.exists(PATH) else {}
    J = json.load(open(os.path.join(OUT, 'r58_theory.json')))['J_inf_bits']
    if sys.argv[1] == 'data':
        for name, (L, W) in WINS.items():
            if name in res and 'J_data' in res[name]:
                continue
            j, n, cells = data(L, W)
            res.setdefault(name, {}).update({'J_data': j, 'triples': n, 'cells': cells, 'mean_lambda': math.log(L + W / 2), 'J_inf': J})
            json.dump(res, open(PATH, 'w'), indent=1); print(name, res[name], flush=True)
    else:
        for name in WINS:
            lam = res[name]['mean_lambda'] if name in res else math.log(WINS[name][0] + WINS[name][1] / 2)
            jl = law(lam, int(12 * lam) // 2 * 2)
            res.setdefault(name, {}).update({'J_law_discrete': jl})
            json.dump(res, open(PATH, 'w'), indent=1); print(name, res[name], flush=True)
