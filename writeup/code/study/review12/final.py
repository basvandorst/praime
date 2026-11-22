"""Order-2 law by Gauss-Legendre (share at each node) and order 3 + overlap with delta(lam) from a cubic in 1/lam through
the five computed values (o3_*.json); leave-one-out check of that fit at lam = 34 and 40.5."""
import math, json
import numpy as np
from law import Gap, primes_upto

data = json.load(open('/home/bas/praime/writeup/data/study/oes_gap_counts_4e18.json'))['N']
f = json.load(open('/home/bas/praime/writeup/data/study/r53f.json'))
o3 = {**json.load(open('o3_72_300.json')), **json.load(open('o3_700.json'))}
P = primes_upto(10_000_100); P = P[P <= 10_000_019]; g0 = np.diff(P)
A, B = math.log(1e7), math.log(4e18)
edges = [A, 26, 32, 36, 38, 39.5, 40.5, 41.3, 42.0, B]
x, wt = np.polynomial.legendre.leggauss(10)
out = {}
for g in (72, 300, 700):
    G = Gap(g); lams = np.array(o3[str(g)]['lams']); d = np.array(o3[str(g)]['delta_own'])
    cf = np.polyfit(1 / lams, d, 3)
    loo = {}
    for k in (1, 3):
        m = np.ones(5, bool); m[k] = False
        c2 = np.polyfit(1 / lams[m], d[m], 3)
        loo[float(lams[k])] = float(np.polyval(c2, 1 / lams[k]) - d[k])
    I2 = I3 = 0.0
    for lo, hi in zip(edges[:-1], edges[1:]):
        for xi, wi in zip(x, wt):
            u = 0.5 * (hi - lo) * xi + 0.5 * (hi + lo)
            l2 = G.lnshare2(u)
            dl = np.polyval(cf, 1 / max(u, 30.0)) if u >= 30 else np.polyval(cf, 1 / 30.0)
            I2 += 0.5 * (hi - lo) * wi * math.exp(l2 + u) / u
            I3 += 0.5 * (hi - lo) * wi * math.exp(l2 + dl + u) / u
    base = int(np.sum(g0 == g)); D = data[str(g)]
    law2, law3 = I2 + base, I3 + base
    out[g] = {'rel_o2': (D - law2) / D, 'rel_o3': (D - law3) / D, 'step': (law3 - law2) / D,
              'r53f': {k: f[str(g)][k] for k in ('rel_o2', 'rel_o3', 'step')}, 'loo_fit_err': loo}
    print(g, {k: ('%+.4f %%' % (100 * v) if isinstance(v, float) else v) for k, v in out[g].items()}, flush=True)
json.dump(out, open('final.json', 'w'), indent=1)
