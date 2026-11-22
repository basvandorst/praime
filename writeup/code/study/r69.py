"""Round 69 (preregister/round69.txt): conditional information of the gap before last in [1e12, 1e12 + 1e10) against the
discrete leading-order law, with a conditional shuffle control (g3 permuted among triples with the same g2, which
keeps both pair tables and has true value 0) on the first 3e8 numbers, scaled by sample size. usage: r69.py"""
import os, json, math
import numpy as np
import common as C
from r58 import cmi_from_triples, law, SH

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
L0, L1 = 10 ** 12, 10 ** 12 + 10 ** 10


def triples(lo, hi, keep=False):
    acc = {}; prev = None; gs = []
    for s in range(lo, hi, 10 ** 8):
        P = np.nonzero(C.primes_mask(s, 10 ** 8))[0].astype(np.int64) + s
        if prev is not None:
            P = np.concatenate([prev, P])
        g = np.diff(P)
        if keep:
            gs.append(g if prev is None else g[2:])          # the first two gaps repeat the previous segment's last two
        k = (g[:-2] * SH + g[1:-1]) * SH + g[2:]
        v, c = np.unique(k, return_counts=True)
        for a, b in zip(v.tolist(), c.tolist()):
            acc[a] = acc.get(a, 0) + b
        prev = P[-3:]
    keys = np.array(list(acc.keys()), np.int64); cnt = np.array(list(acc.values()), float)
    return keys, cnt, (np.concatenate(gs) if keep else None)


if __name__ == '__main__':
    lam = float(np.mean(np.log(np.linspace(L0, L1, 1001))))
    J = json.load(open(os.path.join(OUT, 'r58_theory.json')))['J_inf_bits']
    keys, cnt, _ = triples(L0, L1)
    Jd = cmi_from_triples(keys, cnt); N = cnt.sum()
    # shuffle control on the first 3e8 numbers
    ks, cs, g = triples(L0, L0 + 3 * 10 ** 8, keep=True)
    del ks, cs
    g = g.astype(np.int32); g1, g2, g3 = g[:-2], g[1:-1], g[2:].copy()
    rng = np.random.default_rng(69)
    order = np.argsort(g2, kind='stable'); g2s = g2[order]
    bounds = np.flatnonzero(np.diff(g2s)) + 1
    g3o = g3[order]
    for blk in np.split(np.arange(len(g2s)), bounds):
        g3o[blk] = rng.permutation(g3o[blk])
    g3s = np.empty_like(g3); g3s[order] = g3o
    k = (g1.astype(np.int64) * SH + g2) * SH + g3s
    v, c = np.unique(k, return_counts=True)
    bias_sub = cmi_from_triples(v, c.astype(float)); Nsub = float(len(k))
    bias = bias_sub * Nsub / N
    Jl = law(lam, int(20 * lam) // 2 * 2)
    out = {'lambda': lam, 'triples': int(N), 'cells': len(cnt), 'J_data': Jd, 'shuffle_bias_sub': bias_sub, 'n_sub': int(Nsub),
           'bias_scaled': bias, 'J_data_corrected': Jd - bias, 'J_law_discrete': Jl, 'J_inf': J,
           'corrected_minus_law': Jd - bias - Jl, 'lambda_times_corrected_minus_Jinf': lam * (Jd - bias - J)}
    print({a: (round(b, 5) if isinstance(b, float) else b) for a, b in out.items()})
    json.dump(out, open(os.path.join(OUT, 'r69.json'), 'w'), indent=1)
