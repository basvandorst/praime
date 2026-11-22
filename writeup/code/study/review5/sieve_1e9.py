"""Own sieve of [1e9, 3e9): gap histogram on [1e9, 2e9) (round 30) and twin counts per 1e6 interval (round 31 addendum 2)."""
import numpy as np, json, math
from lib import seg_mask

hist = np.zeros(2000, np.int64); last = None; nprimes = 0
tw6 = []
for s in range(200):
    lo = 10 ** 9 + s * 10 ** 7
    M = seg_mask(lo, 10 ** 7 + 2)
    tw = (M[:-2] & M[2:])[:10 ** 7]
    tw6 += tw.reshape(10, 10 ** 6).sum(1).tolist()
    if lo < 2 * 10 ** 9:
        P = np.nonzero(M[:10 ** 7])[0].astype(np.int64) + lo
        nprimes += len(P)
        if last is not None:
            P = np.concatenate(([last], P))
        hist += np.bincount(np.minimum(np.diff(P), 1999), minlength=2000)
        last = int(P[-1])
    if s % 20 == 0:
        print(s, flush=True)
np.save('/home/bas/scratch/study-review5/hist_1e9_2e9.npy', hist)
np.save('/home/bas/scratch/study-review5/tw6_1e9_3e9.npy', np.array(tw6))
print('primes in [1e9,2e9):', nprimes, 'gaps:', int(hist.sum()))
for a, b in [(120, 178), (180, 238), (240, 298), (300, 358), (360, 454)]:
    print(a, b, int(hist[a:b + 1].sum()))
print('max gap', int(np.nonzero(hist)[0].max()))
