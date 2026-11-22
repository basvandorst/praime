"""Own MI of consecutive gaps in a window: plug-in, Miller-Madow, shuffle control (true MI 0), half-sample extrapolation."""
import sys, math, json
import numpy as np
sys.path.insert(0, '/home/bas/scratch/study-review13')
from sieve import base_primes, odd_primes_in

def H(c):
    c = c[c > 0]; n = c.sum(); q = c / n
    return float(-np.sum(q * np.log2(q))), len(c), n

def mi(g1, g2, mm=True):
    key = g1.astype(np.int64) * 4096 + g2
    _, cj = np.unique(key, return_counts=True); _, c1 = np.unique(g1, return_counts=True); _, c2 = np.unique(g2, return_counts=True)
    (hj, kj, n), (h1, k1, _), (h2, k2, _) = H(cj), H(c1), H(c2)
    I = h1 + h2 - hj
    corr = (kj - k1 - k2 + 1) / (2 * n * math.log(2))
    return I, I - corr, kj

L, W = int(float(sys.argv[1])), int(float(sys.argv[2]))
bp = base_primes(int(math.isqrt(L + W)) + 2)
P = np.concatenate([odd_primes_in(s, min(s + 10 ** 7, L + W), bp) for s in range(L, L + W, 10 ** 7)])
g = np.diff(P); g1, g2 = g[:-1], g[1:]
print('primes', len(P), 'pairs', len(g1), 'mean lam %.3f' % np.mean(np.log(P.astype(float))))
pi, mmI, kj = mi(g1, g2); print('plug-in %.4f  MM %.4f  cells %d  MM correction %.4f' % (pi, mmI, kj, pi - mmI))
rng = np.random.default_rng(1)
sh = [mi(g1, rng.permutation(g2))[1] for _ in range(5)]
print('shuffled (true 0): MM mean %.5f sd %.5f' % (np.mean(sh), np.std(sh)))
h = len(g1) // 2
a = mi(g1[:h], g2[:h])[1]; b = mi(g1[h:], g2[h:])[1]
print('halves MM %.4f %.4f; 2*full - mean(halves) %.4f' % (a, b, 2 * mmI - (a + b) / 2))
# block-bootstrap sd: 20 blocks
blocks = np.array_split(np.arange(len(g1)), 20)
vals = [mi(np.concatenate([g1[blocks[j]] for j in range(20) if j != i]), np.concatenate([g2[blocks[j]] for j in range(20) if j != i]))[1] for i in range(20)]
print('jackknife sd %.4f' % (math.sqrt(19 / 20 * np.sum((np.array(vals) - np.mean(vals)) ** 2))))
