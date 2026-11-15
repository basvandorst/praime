"""Extend trial factoring of M_p = 2^p - 1 from depth K0 to depth K (factors q = 2kp + 1, k <= K) for the primes p < X
that have no factor with k <= K0 in ../kmin_<X>_<K0>.npy.
Per p, vectorised over k in chunks: keep q = +-1 mod 8, sieve out q with a prime factor l <= 3000
(q = 0 mod l  <=>  k = -(2p)^-1 mod l; q > l always since 2(K0+1)p + 1 > 3000), then 2^p mod q by left-to-right
binary powering with a float-assisted uint64 mulmod (exact for q < 2^50; q < 2K X + 1 here).
Output: kmin_<X>_<K>.npy = [P, kmin] (kmin = smallest k <= K with q | M_p, 0 = none), plus the count of candidate q
that reached the powering step (the TF work) per p in tfwork_<X>_<K>.npy."""
import sys, time, numpy as np

X, K0, K = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
import os
P, kmin = np.load(f'kmin_{X}_{K0}.npy' if os.path.exists(f'kmin_{X}_{K0}.npy') else f'../kmin_{X}_{K0}.npy')
work0 = np.load(f'tfwork_{X}_{K0}.npy') if os.path.exists(f'tfwork_{X}_{K0}.npy') else None
assert 2 * K * X + 1 < 2 ** 50
SP = [l for l in range(3, 3001) if all(l % d for d in range(2, int(l ** 0.5) + 1))]
assert 2 * (K0 + 1) * 3 + 1 > 3000
C = 1 << 20

def pow2_mod(p, q):
    """2^p mod q, elementwise over uint64 array q (q < 2^50), scalar exponent p."""
    qi = q.view(np.int64); qf = q.astype(np.float64)
    r = np.ones(len(q), dtype=np.uint64)
    for b in bin(p)[2:]:
        h = np.floor(r.astype(np.float64) * r.astype(np.float64) / qf).astype(np.uint64)
        t = (r * r - h * q).view(np.int64)            # wraps mod 2^64; true value in [-q, 2q)
        t = np.where(t < 0, t + qi, t); t = np.where(t >= qi, t - qi, t)
        if b == '1':
            t = t << 1; t = np.where(t >= qi, t - qi, t)
        r = t.view(np.uint64)
    return r

work = np.zeros(len(P), dtype=np.int64) if work0 is None else work0.copy()   # TF candidates counted from the first file's K0
t0 = time.time()
todo = np.flatnonzero(kmin == 0)
for n, i in enumerate(todo):
    p = int(P[i])
    inv = {l: (-pow(2 * p, -1, l)) % l for l in SP if l != p}
    for s in range(K0 + 1, K + 1, C):
        e = min(s + C, K + 1)
        k = np.arange(s, e, dtype=np.uint64)
        q = np.uint64(2 * p) * k + np.uint64(1)
        m8 = q & np.uint64(7)
        keep = (m8 == 1) | (m8 == 7)
        for l, k0 in inv.items():
            keep[(k0 - s) % l::l] = False
        kk, qq = k[keep], q[keep]
        work[i] += len(qq)
        hit = np.flatnonzero(pow2_mod(p, qq) == 1)
        if p < 64: hit = np.array([h for h in hit if int(qq[h]) < (1 << p) - 1], dtype=np.int64)   # proper divisors only
        if len(hit):
            kmin[i] = int(kk[hit[0]]); break
    if n % 500 == 0: print(f'{n}/{len(todo)} p={p} {time.time() - t0:.0f} s', flush=True)
print(f'X {X} K {K}: {len(P)} odd primes, {np.sum(kmin == 0)} without a factor 2kp+1 (k <= {K}), {time.time() - t0:.0f} s')
np.save(f'kmin_{X}_{K}.npy', np.stack([P, kmin]))
np.save(f'tfwork_{X}_{K}.npy', work)
