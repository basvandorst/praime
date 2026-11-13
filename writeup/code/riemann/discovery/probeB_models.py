"""Waiting times for the first prime along a progression a, a + q, a + 2q, ... (q = 30,000,001), against models.
Needs probeB_30000001_j.npy from probeB.py. Model with sieve level z: a member with no prime factor <= z is prime with
probability prod_{p<=z} p/(p-1) / ln m, otherwise 0 (z = 2 is independence apart from parity). Then the champions:
how many of their first members escape the primes below 1000, against random classes."""
import numpy as np, math
q = 30000001
J = np.load('probeB_30000001_j.npy').astype(np.int64)
samp = np.random.default_rng(1).choice(np.arange(1, q, dtype=np.int64), 200000, replace=False)
def primes_upto(z): return [p for p in range(2, z + 1) if all(p % r for r in range(2, int(p ** 0.5) + 1))]
js = [10, 50, 100, 150, 200]
data = {j: float(np.mean(J[1:] >= j)) for j in js}
for z in (2, 13, 31, 101, 307):
    P = primes_upto(z); boost = np.prod([p / (p - 1) for p in P]); logs = np.zeros(len(samp)); res = {}
    for i in range(max(js)):
        m = samp + i * q; ok = np.ones(len(m), dtype=bool)
        for p in P: ok &= (m % p != 0)
        logs += np.log1p(-np.where(ok, np.minimum(boost / np.log(m.astype(float)), 0.999), 0.0))
        if i + 1 in js: res[i + 1] = float(np.mean(np.exp(logs)))
    print(f'z = {z:3d}: data / model ' + ', '.join(f'j={j}: {data[j] / res[j]:.3f}' for j in js))
P = primes_upto(1000)[1:]
def survivors(a, n):
    m = a + q * np.arange(n, dtype=np.int64); ok = (m % 2 == 1)
    for p in P: ok &= (m % p != 0)
    return int(ok.sum())
top = np.argsort(J)[::-1][:20]
s_top = [survivors(int(a), 313) for a in top]
s_rnd = np.array([survivors(int(a), 313) for a in np.random.default_rng(5).choice(np.arange(1, q), 2000, replace=False)])
print(f'first 313 members free of primes below 1000: 20 latest classes {s_top} (mean {np.mean(s_top):.1f}); random classes mean {s_rnd.mean():.1f}, sd {s_rnd.std():.1f}')
