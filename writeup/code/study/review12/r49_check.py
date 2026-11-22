"""Round 49 spot check: law for g = 100, 140, 200 in [1e7, 1e8) by Gauss-Legendre (no interpolation), order 2, and
order 3 + overlap with delta at 5 lam; against r49.json (3-point linear interpolation of ln share). Counts recomputed."""
import math, json
import numpy as np, mpmath
from law import Gap, primes_upto

r49 = json.load(open('/home/bas/praime/writeup/data/study/r49.json'))
L0, L1 = 1e7, 1e8
A, B = math.log(L0), math.log(L1)
x, wt = np.polynomial.legendre.leggauss(16)
norm = r49['norm']
# counts by segments
base = primes_upto(10 ** 4 + 10)
def seg(lo, hi):
    m = np.ones(hi - lo, bool)
    for p in base:
        if p * p >= hi: break
        s = max(p * p, ((lo + p - 1) // p) * p); m[s - lo::p] = False
    return np.nonzero(m)[0] + lo
cnt = {}; last = None
for lo in range(10 ** 7, 10 ** 8, 10 ** 7):
    P = seg(lo, lo + 10 ** 7)
    if last is not None: P = np.concatenate([[last], P])
    last = int(P[-1]); g = np.diff(P)
    for v in (100, 140, 200): cnt[v] = cnt.get(v, 0) + int(np.sum(g == v))
for g in (100, 140, 200):
    G = Gap(g); lams = [16.12, 16.7, 17.2, 17.8, 18.42]; d = G.delta3(lams); cf = np.polyfit(1 / np.array(lams), d, 3)
    I2 = I3 = 0.0
    for xi, wi in zip(x, wt):
        u = 0.5 * (B - A) * xi + 0.5 * (B + A); f = 0.5 * (B - A) * wi * math.exp(G.lnshare2(u) + u) / u
        I2 += f; I3 += f * math.exp(np.polyval(cf, 1 / u))
    pg = r49['per_gap'][str(g)]
    print(g, 'count', cnt[g], pg['count'], 'law_o2 own %.2f r49 %.2f (%+.2f %%)' % (I2 * norm, pg['law_o2'], 100 * (I2 * norm / pg['law_o2'] - 1)),
          'law_o3 own %.2f r49 %.2f (%+.2f %%)' % (I3 * norm, pg['law_o3'], 100 * (I3 * norm / pg['law_o3'] - 1)), flush=True)
