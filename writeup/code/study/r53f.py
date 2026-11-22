"""Round 53, corrected numerics: the coarse lam grid of r53.py (steps of 2, linear interpolation of ln share) biased
the integral low by an amount growing with g (ln share contains -g/lam, curvature 2g/lam^3); r53e.py shows it. Here
order 2 on a fine grid (step 1 below 30, 0.25 from 30 to ln 4e18) plus the order-3 + overlap correction (lam = 30, 38,
42.84 as in r53.py; it varies slowly), for g = 40..70 (where the exact law fits the counts to about 1e-8, round 23)
and the 20 gaps of round 53. usage: r53f.py"""
import os, json, math, time
import numpy as np
import common as C
import r47
from r51 import o3_delta
from r53 import GS, L3

data = json.load(open(os.path.join(r47.OUT, 'oes_gap_counts_4e18.json')))['N']
old = json.load(open(os.path.join(r47.OUT, 'r53.json')))
M = C.primes_mask(2, 10 ** 7 - 2); P = np.append(np.nonzero(M)[0] + 2, 10000019); g0 = np.diff(P)
u = np.linspace(math.log(1e7), math.log(4e18), 40001)
LF = np.concatenate([np.arange(16.0, 30.0, 1.0), np.arange(30.0, 42.84, 0.25), [math.log(4e18)]])
path = os.path.join(r47.OUT, 'r53f.json')
out = json.load(open(path)) if os.path.exists(path) else {}
t0 = time.time()
for g in [40, 50, 60, 70] + GS:
    if str(g) in out or time.time() - t0 > 450:
        continue
    st = r47.setup(g); base = int(np.sum(g0 == g)); d = data[str(g)]
    ls2 = np.interp(u, LF, [math.log(r47.share(st, lam, 2)) for lam in LF])
    d3 = old['o3'][str(g)] if str(g) in old['o3'] else o3_delta(g, L3)
    ls3 = ls2 + np.interp(u, L3, d3)
    law2 = float(np.trapezoid(np.exp(ls2) * np.exp(u) / u, u)) + base
    law3 = float(np.trapezoid(np.exp(ls3) * np.exp(u) / u, u)) + base
    out[str(g)] = {'data': d, 'rel_o2': (d - law2) / d, 'rel_o3': (d - law3) / d, 'step': (law3 - law2) / d, 'sd_rel': 1 / math.sqrt(d),
                   'pass': bool(abs((d - law3) / d) <= max(0.001, 2 * abs((law3 - law2) / d)))}
    json.dump(out, open(path, 'w'), indent=1)
    print(g, {k: ('%+.4f %%' % (100 * v) if isinstance(v, float) and k != 'sd_rel' else v) for k, v in out[str(g)].items() if k != 'data'}, flush=True)
print('done', len(out), flush=True)
