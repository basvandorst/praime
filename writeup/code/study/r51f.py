"""Round 51 with the fine lam grid (r53e.py showed that linear interpolation of ln share on steps of 1-2 biases the law
low by an amount growing with g): order 2 on steps of 0.25 from 16 to 27, plus round 51's order-3 + overlap
corrections. usage: r51f.py (resumable)"""
import os, json, math, time
import numpy as np
from scipy.stats import kstest
import r47
from r51 import GS, P

path = os.path.join(r47.OUT, 'r51f.json')
r51 = json.load(open(os.path.join(r47.OUT, 'r51.json')))
LF = np.arange(16.0, 27.01, 0.25)
res = json.load(open(path)) if os.path.exists(path) else {'o2': {}}
t0 = time.time()
for g in GS:
    if str(g) in res['o2'] or time.time() - t0 > 450:
        continue
    st = r47.setup(g); res['o2'][str(g)] = [math.log(r47.share(st, lam, 2)) for lam in LF]
    json.dump(res, open(path, 'w'))
done = [g for g in GS if str(g) in res['o2']]
print('done', len(done), flush=True)
if len(done) == len(GS):
    out = {}
    for name, use3 in (('order2', False), ('order3_overlap', True)):
        U = []
        for g in GS:
            u = np.linspace(math.log(1e7), math.log(P[g]), 8001)
            ls = np.interp(u, LF, res['o2'][str(g)])
            if use3:
                o = r51['o3'][str(g)]; ls = ls + np.interp(u, o['lams'], o['delta'])
            U.append(float(np.trapezoid(np.exp(ls) / u * np.exp(u), u)))
        U = np.array(U); garr = np.array(GS); ks = kstest(U, 'expon')
        out[name] = {'mean': float(U.mean()), 'z_null': float((U.mean() - 1) * math.sqrt(len(U))), 'ks_p': float(ks.pvalue),
                     'mean_180_300': float(U[garr <= 300].mean()), 'mean_302_454': float(U[garr > 300].mean())}
        print(name, {k: round(v, 4) for k, v in out[name].items()}, flush=True)
    res['summary'] = out
    json.dump(res, open(path, 'w'), indent=1)
