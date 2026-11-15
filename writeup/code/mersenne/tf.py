"""Experiment 1: hit rate of Mersenne exponents among primes p < X after trial factoring 2^p - 1.
A factor of M_p = 2^p - 1 has the form q = 2kp + 1 with q = +-1 mod 8. For every prime p < X find the smallest
k <= K with q | M_p (q < M_p), vectorised over p: q < 2^31, so q^2 < 2^62 fits in uint64.
Output: kmin per prime (0 = none up to K) in kmin_<X>_<K>.npy."""
import sys, numpy as np, time
from sympy import primerange
X, K = int(sys.argv[1]), int(sys.argv[2])
P = np.array(list(primerange(3, X)), dtype=np.uint64)
assert 2 * K * int(P.max()) + 1 < 2 ** 31
kmin = np.zeros(len(P), dtype=np.int64)
nb = int(P.max()).bit_length()
t0 = time.time()
for k in range(1, K + 1):
    open_ = kmin == 0
    if not open_.any(): break
    p = P[open_]; q = 2 * np.uint64(k) * p + np.uint64(1)
    ok = (q % np.uint64(8) == 1) | (q % np.uint64(8) == 7)
    if not ok.any(): continue
    p, q = p[ok], q[ok]
    r = np.ones(len(q), dtype=np.uint64)
    for b in range(nb - 1, -1, -1):          # left-to-right square and multiply, exponent p (per element)
        r = (r * r) % q
        bit = (p >> np.uint64(b)) & np.uint64(1)
        r = np.where(bit == 1, (r * np.uint64(2)) % q, r)
    hit = r == 1
    # q must be a proper divisor: q < M_p (only p = 3, k = 1 gives q = 7 = M_3)
    small = p < 64
    proper = np.array([int(qq) < (1 << int(pp)) - 1 for pp, qq in zip(p, q)]) if small.any() else np.ones(len(p), bool)
    hit &= np.where(small, proper, True)
    idx = np.flatnonzero(open_)[ok][hit]
    kmin[idx] = k
print(f'X {X} K {K}: {len(P)} odd primes, {np.sum(kmin == 0)} without a factor 2kp+1 (k <= {K}), {time.time() - t0:.0f} s')
np.save(f'kmin_{X}_{K}.npy', np.stack([P.astype(np.int64), kmin]))
