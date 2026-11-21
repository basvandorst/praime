"""Round 3 statistics of the true primes, windows train, A, B and C, whole and in 10 blocks.
Writes data/study/results_r3.jsonl: the round-2 statistics (probes.sequential) plus tails, CPAP-4, entropies and the
19x19 cell counts of consecutive gap pairs (h1, h2 <= 38) with the total number of pairs."""
import json, os, math
import numpy as np
import common as C
import probes as PR
from theory_r3 import WIN, tail_thresholds

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')


def H_bits(c):
    p = c / c.sum(); p = p[p > 0]
    return float(-np.sum(p * np.log2(p)))


def stats(M, L, W, win):
    d = PR.sequential(M, L, W)
    P = np.nonzero(M)[0].astype(np.int64) + L
    g = np.diff(P)
    n = len(g)
    for th in tail_thresholds(win):
        d['T1.tail%d' % th] = float(np.mean(g >= th))
    d['T3.cpap4'] = float(np.mean((g[2:] == g[1:-1]) & (g[1:-1] == g[:-2]))) * (n - 2) / n
    c = np.minimum(g // 2, 20)
    d['T4.H_class'] = H_bits(np.bincount(c, minlength=21)[1:].astype(float))
    pair = np.bincount((c[:-1] - 1) * 20 + (c[1:] - 1), minlength=400).astype(float)
    d['T4.H_cond'] = H_bits(pair) - H_bits(np.bincount(c[:-1], minlength=21)[1:].astype(float))
    small = (g[:-1] <= 38) & (g[1:] <= 38)
    cells = np.bincount((g[:-1][small] // 2 - 1) * 19 + (g[1:][small] // 2 - 1), minlength=361)
    d['_cells'] = cells.tolist()
    d['_npairs'] = int(n - 1)
    return d


def main():
    path = os.path.join(OUT, 'results_r3.jsonl')
    with open(path, 'w') as f:
        for win, (L, W) in WIN.items():
            M = C.primes_mask(L, W)
            for b in range(-1, 10):
                Lb, Wb = (L, W) if b < 0 else (L + b * (W // 10), W // 10)
                Mb = M[Lb - L:Lb - L + Wb]
                f.write(json.dumps({'round': 'r3', 'window': win, 'gen': 'primes', 'block': b, 'L': Lb, 'W': Wb,
                                    'members': int(Mb.sum()), 'stats': stats(Mb, Lb, Wb, win)}) + '\n')
            print(win, flush=True)


if __name__ == '__main__':
    main()
