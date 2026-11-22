"""Round 52 (preregister/round52.txt): entropy of the next prime gap against log2(e lam) - K/ln 2 (leading order of
the Hardy-Littlewood law). The formula's entropy for a window is that of the mixture over the window, P_mix(g) =
sum_t w_t (S(0,g)/lam_t) e^(-g/lam_t) / sum_t w_t with w_t the prime density, so window width is handled exactly.
usage: r52.py"""
import os, json, math
import numpy as np
import common as C
import theory_r3 as R3

PR = C.small_primes(20_000_000).astype(np.float64)
P3 = PR[PR > 2]
C2 = math.exp(float(np.sum(np.log1p(-1 / (P3 - 1) ** 2))))
K = math.log(2 * C2) + float(np.sum(np.log((P3 - 1) / (P3 - 2)) / (P3 - 1))) + 1 / (2e7 * math.log(2e7))
WINS = {'train': (10 ** 7, 9 * 10 ** 7), '1e9-2e9': (10 ** 9, 10 ** 9), 'A': (10 ** 10, 5 * 10 ** 7),
        'B': (10 ** 12, 5 * 10 ** 7), 'C': (10 ** 15, 2 * 10 ** 7)}


def S_even(G):
    k = np.arange(1, G // 2 + 1); S = np.full(len(k), 2 * C2)
    for p in C.small_primes(G // 2 + 1)[1:]:
        p = int(p); S[p - 1::p] *= (p - 1) / (p - 2)
    return S                                                     # S(0, 2k), k = 1..G/2


def H_formula(L, W, G=4000):
    S = S_even(G); g = 2 * np.arange(1, G // 2 + 1)
    t = np.linspace(L, L + W, 201); lam = np.log(t); w = 1 / lam
    Pm = np.zeros(len(g))
    for li, wi in zip(lam, w):
        Pm += wi * S / li * np.exp(-g / li)
    Pm /= Pm.sum()
    return float(-np.sum(Pm * np.log2(Pm))), float(np.average(lam, weights=w))


def H_data(L, W):
    cnt = {}; last = None; seg = 10 ** 7
    for s in range(L, L + W, seg):
        P = np.nonzero(C.primes_mask(s, min(seg, L + W - s)))[0].astype(np.int64) + s
        if last is not None:
            P = np.concatenate([[last], P])
        last = int(P[-1])
        v, c = np.unique(np.diff(P), return_counts=True)
        for a, b in zip(v.tolist(), c.tolist()):
            cnt[a] = cnt.get(a, 0) + b
    c = np.array(list(cnt.values()), dtype=float); n = c.sum(); p = c / n
    return float(-np.sum(p * np.log2(p)) + (len(c) - 1) / (2 * n * math.log(2))), int(n)


if __name__ == '__main__':
    out = {'K_nats': K, 'K_minus_ln2_bits': (K - math.log(2)) / math.log(2)}
    for name, (L, W) in WINS.items():
        hd, n = H_data(L, W); hf, lam = H_formula(L, W)
        asym = (math.log(math.e * lam) - K) / math.log(2)
        out[name] = {'gaps': n, 'mean_lambda': lam, 'H_data': hd, 'H_formula_mixture': hf, 'H_asym_at_mean_lambda': asym,
                     'diff': hd - hf, 'lambda_times_diff': lam * (hd - hf)}
        print(name, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in out[name].items()}, flush=True)
    json.dump(out, open(os.path.join(R3.OUT, 'r52.json'), 'w'), indent=1)
