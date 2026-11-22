"""Round 53 (preregister/round53.txt): the conditioned expansion against the counts of gaps below 4e18.
usage: r53.py o2 | o3 | compare   (o2/o3 resumable)"""
import os, sys, json, math, time
import numpy as np
import common as C
import r47
from r51 import o3_delta

OUT = r47.OUT
PATH = os.path.join(OUT, 'r53.json')
GS = [72, 100, 150, 200, 250, 300, 350, 400, 454, 460, 482, 500, 520, 538, 556, 580, 610, 640, 670, 700]
L2 = np.concatenate([np.arange(16.0, 42.0, 2.0), [42.0, 42.84]])
L3 = np.array([30.0, 38.0, 42.84])

if __name__ == '__main__':
    res = json.load(open(PATH)) if os.path.exists(PATH) else {'o2': {}, 'o3': {}}
    t0 = time.time(); mode = sys.argv[1]
    if mode in ('o2', 'o3'):
        for g in GS:
            if str(g) in res[mode] or time.time() - t0 > 380:
                continue
            if mode == 'o2':
                st = r47.setup(g); res['o2'][str(g)] = [math.log(r47.share(st, lam, 2)) for lam in L2]
            else:
                res['o3'][str(g)] = o3_delta(g, L3)
            json.dump(res, open(PATH, 'w'), indent=1)
            print(mode, g, round(time.time() - t0), res[mode][str(g)][-1], flush=True)
        print(mode, 'done', sum(str(g) in res[mode] for g in GS), 'of', len(GS), flush=True)
    else:
        data = json.load(open(os.path.join(OUT, 'oes_gap_counts_4e18.json')))['N']
        old = json.load(open(os.path.join(OUT, 'r29.json')))
        M = C.primes_mask(2, 10 ** 7 - 2); P = np.append(np.nonzero(M)[0] + 2, 10000019); g0 = np.diff(P)
        u = np.linspace(math.log(1e7), math.log(4e18), 20001)
        out = {}
        for g in GS:
            if str(g) not in res['o2']:
                continue
            base = int(np.sum(g0 == g)); d = data[str(g)]
            ls2 = np.interp(u, L2, res['o2'][str(g)])
            law2 = float(np.trapezoid(np.exp(ls2) * np.exp(u) / u, u)) + base
            row = {'data': d, 'law_o2': law2, 'rel_o2': (d - law2) / d, 'sd_rel': 1 / math.sqrt(d)}
            if str(g) in res['o3']:
                ls3 = ls2 + np.interp(u, L3, res['o3'][str(g)])
                law3 = float(np.trapezoid(np.exp(ls3) * np.exp(u) / u, u)) + base
                row.update({'law_o3': law3, 'rel_o3': (d - law3) / d, 'step': (law3 - law2) / d})
            if str(g) in old:
                row['rel_old'] = old[str(g)]['rel']
            out[str(g)] = row
            print(g, {k: ('%.3e' % v if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
        res['compare'] = out
        json.dump(res, open(PATH, 'w'), indent=1)
