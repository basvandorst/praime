"""Round 53 addendum: the gaps 200, 250, 300 with y = 17 (92160 classes) instead of 13, to see whether the ~0.1 %
remainder below 4e18 is truncation of the expansion. Run with R47_Y=17. usage: r53b.py"""
import os, json, math
import numpy as np
import common as C
import r47
from r51 import o3_delta
from r53 import L2, L3

assert r47.Y == 17
data = json.load(open(os.path.join(r47.OUT, 'oes_gap_counts_4e18.json')))['N']
M = C.primes_mask(2, 10 ** 7 - 2); P = np.append(np.nonzero(M)[0] + 2, 10000019); g0 = np.diff(P)
u = np.linspace(math.log(1e7), math.log(4e18), 20001)
out = {}
for g in (200, 250, 300):
    st = r47.setup(g)
    l2 = [math.log(r47.share(st, lam, 2)) for lam in L2]
    d3 = o3_delta(g, L3)
    base = int(np.sum(g0 == g)); d = data[str(g)]
    ls2 = np.interp(u, L2, l2); ls3 = ls2 + np.interp(u, L3, d3)
    law2 = float(np.trapezoid(np.exp(ls2) * np.exp(u) / u, u)) + base
    law3 = float(np.trapezoid(np.exp(ls3) * np.exp(u) / u, u)) + base
    out[str(g)] = {'rel_o2_y17': (d - law2) / d, 'rel_o3_y17': (d - law3) / d, 'step_y17': (law3 - law2) / d}
    print(g, {k: '%+.4f %%' % (100 * v) for k, v in out[str(g)].items()}, flush=True)
json.dump(out, open(os.path.join(r47.OUT, 'r53b_y17.json'), 'w'), indent=1)
