"""Own segmented sieve: the primes in the validation stretch [1.5e11, 1.5e11 + 1e8 + 1e4) only.
Writes data/val_primes.npy (int64). The extra 1e4 lets p_{n+k} exist for every p_n < 1.5e11 + 1e8."""
import numpy as np, math, os
X0 = 150_000_000_000; SPAN = 100_000_000; EXTRA = 10_000
HI = X0 + SPAN + EXTRA
def small_primes(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]
BASE = small_primes(math.isqrt(HI) + 1)
out = []; SEG = 10_000_000
lo = X0
while lo < HI:
    hi = min(lo + SEG, HI)
    seg = np.ones(hi - lo, bool)
    for p in BASE.tolist():
        st = (-lo) % p
        seg[st::p] = False
    out.append(lo + np.nonzero(seg)[0])
    lo = hi
P = np.concatenate(out).astype(np.int64)
os.makedirs('data', exist_ok=True)
np.save('data/val_primes.npy', P)
print('primes in [%d, %d): %d   first %d last %d' % (X0, HI, len(P), P[0], P[-1]))
print('in [X0, X0+SPAN): %d' % np.searchsorted(P, X0 + SPAN))
