"""Independent check of the trial-factoring data (plain integers, no numpy, no sieve, no mod-8 rule):
smallest k with 2^p = 1 mod (2kp + 1) and 2kp + 1 < M_p, by direct search.
  (a) every odd prime p < 10^4, k <= 10^4;  (b) every odd prime p < 10^5, k <= 10^3;
  (c) 40 random p < 10^5 that survive k <= 10^4, k <= 10^6 (gmpy2 powmod), against kmin_100000_1000000.npy."""
import random, numpy as np, gmpy2
from llind import primes_below

def kfirst(p, K):
    M = (1 << p) - 1
    for k in range(1, K + 1):
        q = 2 * k * p + 1
        if q >= M: return 0
        if gmpy2.powmod(2, p, q) == 1: return k
    return 0

for X, Kf, K in [(10000, 'kmin_10000_10000000.npy', 10 ** 4), (100000, 'kmin_100000_1000000.npy', 10 ** 3)]:
    P, km = np.load(Kf)
    d = {int(p): int(k) for p, k in zip(P, km)}
    bad = [p for p in primes_below(X)[1:] if kfirst(p, K) != (d[p] if 0 < d[p] <= K else 0)]
    print(f'(X={X}, k<={K}) all {len(d)} odd primes: mismatches {bad}')
P, km = np.load('kmin_100000_1000000.npy')
pool = [(int(p), int(k)) for p, k in zip(P, km) if k == 0 or k > 10 ** 4]
random.seed(1)
smp = random.sample(pool, 40)
bad = [(p, k) for p, k in smp if kfirst(p, 10 ** 6) != k]
print(f'(c) 40 random survivors of k <= 10^4, searched to k = 10^6: mismatches {bad}; '
      f'{sum(1 for p, k in smp if k)} of them have a factor with 10^4 < k <= 10^6')
