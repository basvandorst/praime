"""Round 5 statistics: the even n in each window with n^2 + 1 prime (gmpy2.is_prime, no sieve), whole and 10 blocks."""
import json, os
import numpy as np
import gmpy2
from theory_landau import FWIN, GS, QS

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')


def stats(ns):
    g = np.diff(ns)
    d = {'count': int(len(ns))}
    for h in GS:
        d['L.gap%d' % h] = float(np.mean(g == h))
    d['L.ge60'] = float(np.mean(g >= 60)); d['L.ge100'] = float(np.mean(g >= 100))
    for q in QS:
        r = ns % q
        d['L.same%d' % q] = float(np.mean(r[1:] == r[:-1]))
    return d


with open(os.path.join(OUT, 'results_landau.jsonl'), 'w') as f:
    for win, (N0, W) in FWIN.items():
        start = N0 + (N0 % 2)
        ns = np.array([n for n in range(start, N0 + W, 2) if gmpy2.is_prime(n * n + 1)], dtype=np.int64)
        for b in range(-1, 10):
            lo, hi = (N0, N0 + W) if b < 0 else (N0 + b * (W // 10), N0 + (b + 1) * (W // 10))
            f.write(json.dumps({'window': win, 'block': b, 'stats': stats(ns[(ns >= lo) & (ns < hi)])}) + '\n')
        print(win, len(ns), flush=True)
