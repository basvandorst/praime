"""Independent recount of twin starts in [1e9, 2e9 + 4000] (own odd-only segmented sieve, nothing from the study).
Saves the sorted twin starts p (p, p+2 both prime) to twins.npy."""
import math, time, sys
import numpy as np

L0, L1 = 10 ** 9, 2 * 10 ** 9 + 4000
t0 = time.time()
lim = math.isqrt(L1) + 2
s = np.ones(lim + 1, bool); s[:2] = False
for i in range(2, math.isqrt(lim) + 1):
    if s[i]:
        s[i * i::i] = False
base = np.nonzero(s)[0][1:]          # odd primes
SEG = 2 * 10 ** 7                    # numbers per segment (even)
tw = []; prev_last = None
for lo in range(L0, L1, SEG):
    hi = min(lo + SEG, L1)
    # odd numbers lo+1, lo+3, ... (lo even)
    n = (hi - lo) // 2
    odd = np.ones(n, bool)           # odd[i] <-> lo + 1 + 2 i
    for p in base:
        p = int(p)
        st = max(p * p, -(-(lo + 1) // p) * p)
        if st % 2 == 0:
            st += p
        if st >= hi:
            continue
        odd[(st - lo - 1) // 2::p] = False
    P = lo + 1 + 2 * np.nonzero(odd)[0].astype(np.int64)
    if prev_last is not None:
        P = np.concatenate([[prev_last], P])
    t = P[:-1][np.diff(P) == 2]
    tw.append(t); prev_last = int(P[-1])
T = np.concatenate(tw)
np.save('/home/bas/scratch/study-review25/twins.npy', T)
print('twins', len(T), 'in [1e9,2e9):', int(np.sum(T < 2 * 10 ** 9)), 'first', T[:3], 'last<2e9', T[T < 2e9][-3:],
      'time', round(time.time() - t0), flush=True)
