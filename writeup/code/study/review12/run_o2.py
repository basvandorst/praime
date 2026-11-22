"""Order-2 conditioned law integrated by Gauss-Legendre in u = ln t (share evaluated at every node, no interpolation),
plus the actual gaps with lower prime < 1e7; compared with Oliveira e Silva's N(g) below 4e18."""
import sys, math, json, time
import numpy as np
from law import Gap, primes_upto

data = json.load(open('/home/bas/praime/writeup/data/study/oes_gap_counts_4e18.json'))['N']
f = json.load(open('/home/bas/praime/writeup/data/study/r53f.json'))
P = primes_upto(10_000_100); P = P[P <= 10_000_019]   # last prime < 1e7 is 9999991, next 10000019
g0 = np.diff(P)
A, B = math.log(1e7), math.log(4e18)
edges = [A, 26, 32, 36, 38, 39.5, 40.5, 41.3, 42.0, B]

def integral(fun, npt):
    x, wt = np.polynomial.legendre.leggauss(npt)
    tot = 0.0
    for lo, hi in zip(edges[:-1], edges[1:]):
        u = 0.5 * (hi - lo) * x + 0.5 * (hi + lo)
        tot += 0.5 * (hi - lo) * sum(wi * math.exp(fun(ui) + ui) / ui for ui, wi in zip(u, wt))
    return tot

out = {}
for g in map(int, sys.argv[1].split(',')):
    t0 = time.time()
    G = Gap(g)
    base = int(np.sum(g0 == g)); d = data[str(g)]
    I6 = integral(G.lnshare2, 6); I10 = integral(G.lnshare2, 10)
    law = I10 + base
    out[g] = {'base': base, 'law_o2': law, 'rel_o2': (d - law) / d, 'quad_6_vs_10': I6 / I10 - 1,
              'r53f_rel_o2': f[str(g)]['rel_o2'], 'secs': time.time() - t0}
    print(g, out[g], flush=True)
json.dump(out, open('/home/bas/scratch/study-review12/o2_%s.json' % sys.argv[1].replace(',', '_'), 'w'), indent=1)
