"""Round 4 statistics of consecutive twin-prime pairs (preregister/round4.txt), whole windows and 10 blocks."""
import json, os
import numpy as np
import common as C
from theory_twins import WIN, W1

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')


def stats(M, L):
    P = np.nonzero(M)[0].astype(np.int64) + L
    T = P[:-1][np.diff(P) == 2]
    T = T[T > 5]
    h = np.diff(T)
    d = {'n_twins': int(len(T)), 'n_gaps': int(len(h))}
    for x in W1:
        d['W1.gap%d' % x] = float(np.mean(h == x))
    small = h <= 300
    d['W2.le300'] = float(np.mean(small))
    for q in (5, 7):
        d['W3.same%d' % q] = float(np.mean((h[small] % q) == 0))
    return d


with open(os.path.join(OUT, 'results_twins.jsonl'), 'w') as f:
    for win, (L, W) in WIN.items():
        M = C.primes_mask(L, W + 400)          # 400 extra numbers so that p + 2 of the last pair is seen
        for b in range(-1, 10):
            Lb, Wb = (L, W) if b < 0 else (L + b * (W // 10), W // 10)
            Mb = M[Lb - L:Lb - L + Wb + 2]
            d = stats(Mb, Lb)
            f.write(json.dumps({'window': win, 'block': b, 'stats': d}) + '\n')
        print(win, d['n_twins'], flush=True)
