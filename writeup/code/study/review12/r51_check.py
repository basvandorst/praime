"""Round 51 spot check: U(g) = expected number of gaps g with lower prime in [1e7, P_first(g)] under the conditioned law,
order 2 by Gauss-Legendre (no interpolation) and order 3 + overlap with delta from (a) linear interpolation on
[18, 22, 26] as r51.py, (b) a cubic in 1/lam through 18, 20, 22, 24, 26 (extrapolated beyond 26 by the same cubic)."""
import sys, math, json
import numpy as np
from law import Gap

P = {int(g): v for g, v in json.load(open('/home/bas/praime/writeup/data/study/oes_first_occurrence.json'))['P'].items()}
r51 = json.load(open('/home/bas/praime/writeup/data/study/r51.json'))
GS = sorted(g for g in P if 180 <= g <= 454)
lams = [18.0, 20.0, 22.0, 24.0, 26.0]
x, wt = np.polynomial.legendre.leggauss(12)
out = {}
for g in map(int, sys.argv[1].split(',')):
    G = Gap(g); d = G.delta3(lams)
    cf = np.polyfit(1 / np.array(lams), d, 3)
    A, B = math.log(1e7), math.log(P[g])
    edges = np.linspace(A, B, 7)
    U2 = Ul = Us = 0.0
    for lo, hi in zip(edges[:-1], edges[1:]):
        for xi, wi in zip(x, wt):
            u = 0.5 * (hi - lo) * xi + 0.5 * (hi + lo); f = 0.5 * (hi - lo) * wi * math.exp(G.lnshare2(u) + u) / u
            U2 += f
            Ul += f * math.exp(np.interp(u, [18.0, 22.0, 26.0], [d[0], d[2], d[4]]))
            Us += f * math.exp(np.polyval(cf, 1 / u))
    i = GS.index(g)
    out[g] = {'lnP': B, 'delta': d, 'r51_delta_at_18_22_26': r51['o3'][str(g)]['delta'] if 'delta' in r51['o3'][str(g)] else r51['o3'][str(g)],
              'U_o2': U2, 'U_o3_linear': Ul, 'U_o3_smooth': Us,
              'r51_first_run_U_o2': r51['summary']['order2']['U'][i], 'r51_first_run_U_o3': r51['summary']['order3_overlap']['U'][i]}
    print(g, json.dumps(out[g]), flush=True)
json.dump(out, open('/home/bas/scratch/study-review12/r51_%s.json' % sys.argv[1].replace(',', '_'), 'w'), indent=1)
