"""Round 18 statistics: the n in each window with n^2 + n + 41 prime (gmpy2.is_prime, no sieve), whole and 10 blocks."""
import json, os
import numpy as np
import gmpy2
from r18 import FWIN, GS, QS

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')


def stats(ns):
    g = np.diff(ns)
    d = {'count': int(len(ns))}
    for h in GS:
        d['L.gap%d' % h] = float(np.mean(g == h))
    d['L.ge15'] = float(np.mean(g >= 15))
    for q in QS:
        r = ns % q
        d['L.same%d' % q] = float(np.mean(r[1:] == r[:-1]))
    return d


with open(os.path.join(OUT, 'results_euler41.jsonl'), 'w') as f:
    for win, (N0, W) in FWIN.items():
        ns = np.array([n for n in range(N0, N0 + W) if gmpy2.is_prime(n * n + n + 41)], dtype=np.int64)
        for b in range(-1, 10):
            lo, hi = (N0, N0 + W) if b < 0 else (N0 + b * (W // 10), N0 + (b + 1) * (W // 10))
            f.write(json.dumps({'window': win, 'block': b, 'stats': stats(ns[(ns >= lo) & (ns < hi)])}) + '\n')
        print(win, len(ns), flush=True)
