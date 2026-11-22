"""Round 56 (preregister/round56.txt): mutual information between consecutive prime gaps against the closed constant
I_inf = sum_{p>=3} I_p (leading order of the Hardy-Littlewood law). usage: r56.py"""
import os, json, math
import numpy as np
import common as C
import theory_r3 as R3

PR = C.small_primes(20_000_000).astype(np.float64); P = PR[PR > 2]
fA = (1 - 2 / P) / (1 - 1 / P) ** 3; rA = (1 - 1 / P) / (1 - 2 / P)
fB = (1 - 3 / P) / (1 - 1 / P) ** 3; rB = np.where(P > 3, (1 - 3 / P) * (1 - 1 / P) / (1 - 2 / P) ** 2, 1.0)
Ip = (P - 1) / P ** 2 * fA * np.log(rA) + np.where(P > 3, (P - 1) * (P - 2) / P ** 2 * fB * np.log(rB), 0.0)
I_INF = float(Ip.sum()) / math.log(2)
WINS = {'train': (10 ** 7, 9 * 10 ** 7), '1e9-2e9': (10 ** 9, 10 ** 9), 'A': (10 ** 10, 5 * 10 ** 7),
        'B': (10 ** 12, 5 * 10 ** 7), 'C': (10 ** 15, 2 * 10 ** 7)}


def mi(L, W):
    pairs = {}; m1 = {}; last2 = None; seg = 10 ** 7; n = 0
    prev = None
    for s in range(L, L + W, seg):
        Pp = np.nonzero(C.primes_mask(s, min(seg, L + W - s)))[0].astype(np.int64) + s
        if prev is not None:
            Pp = np.concatenate([prev, Pp])
        g = np.diff(Pp)
        key = g[:-1] * 4096 + g[1:]
        v, c = np.unique(key, return_counts=True)
        for a, b in zip(v.tolist(), c.tolist()):
            pairs[a] = pairs.get(a, 0) + b
        prev = Pp[-2:]
    c = np.array(list(pairs.values()), float); k = np.array(list(pairs.keys()))
    N = c.sum(); g1 = k // 4096; g2 = k % 4096
    H = lambda cnt: float(-np.sum((cnt / cnt.sum()) * np.log2(cnt / cnt.sum())) + (len(cnt) - 1) / (2 * cnt.sum() * math.log(2)))
    def marg(keys):
        u, inv = np.unique(keys, return_inverse=True); return np.bincount(inv, weights=c)
    Hj = H(c); H1 = H(marg(g1)); H2 = H(marg(g2))
    return H1 + H2 - Hj, int(N), len(c)


if __name__ == '__main__':
    out = {'I_inf_bits': I_INF}
    print('I_inf bits %.6f' % I_INF, flush=True)
    for name, (L, W) in WINS.items():
        I, N, cells = mi(L, W)
        out[name] = {'pairs': N, 'cells': cells, 'I_bits': I, 'diff': I - I_INF, 'mean_lambda': math.log(L + W / 2)}
        print(name, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in out[name].items()}, flush=True)
    json.dump(out, open(os.path.join(R3.OUT, 'r56.json'), 'w'), indent=1)
