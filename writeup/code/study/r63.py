"""Round 63 (preregister/round63.txt): counts of gaps 276..500 in [1e12, 1e12 + 1e10) against the conditioned law.
usage: r54.py count | o2 | o3 | compare   (resumable)"""
import os, sys, json, math, time
import numpy as np
import mpmath
import common as C
import r47
from r51 import o3_delta

OUT = r47.OUT
PATH = os.path.join(OUT, 'r63.json')
L0, L1 = 10 ** 12, 10 ** 12 + 10 ** 10
GS = list(range(276, 501, 2))
LF = np.linspace(math.log(L0), math.log(L1), 9)
L3 = np.array([math.log(L0), math.log(L1)])
BANDS = [(276, 330), (332, 390), (392, 500)]

if __name__ == '__main__':
    res = json.load(open(PATH)) if os.path.exists(PATH) else {'o2': {}, 'o3': {}}
    mode = sys.argv[1]; t0 = time.time()
    if mode == 'count':
        cnt = {}; last = None; ng = 0
        for s in range(L0, L1, 10 ** 8):
            P = np.nonzero(C.primes_mask(s, 10 ** 8))[0].astype(np.int64) + s
            if last is not None:
                P = np.concatenate([[last], P])
            last = int(P[-1]); g = np.diff(P); ng += len(g)
            v, c = np.unique(g[g >= 200], return_counts=True)
            for a, b in zip(v.tolist(), c.tolist()):
                cnt[str(a)] = cnt.get(str(a), 0) + b
        res['counts'] = cnt; res['gaps'] = ng
        json.dump(res, open(PATH, 'w'), indent=1); print('gaps', ng, 'g>=230:', sum(v for k, v in cnt.items() if int(k) >= 230), flush=True)
    elif mode in ('o2', 'o3'):
        for g in GS:
            if str(g) in res[mode] or time.time() - t0 > 420:
                continue
            if mode == 'o2':
                st = r47.setup(g); res['o2'][str(g)] = [math.log(r47.share(st, lam, 2)) for lam in LF]
            else:
                res['o3'][str(g)] = o3_delta(g, L3)
            json.dump(res, open(PATH, 'w'), indent=1)
        print(mode, 'done', sum(str(g) in res[mode] for g in GS), 'of', len(GS), flush=True)
    else:
        norm = res['gaps'] / float(mpmath.li(L1) - mpmath.li(L0))
        u = np.linspace(math.log(L0), math.log(L1), 4001)
        law = {}
        for g in GS:
            ls = np.interp(u, LF, res['o2'][str(g)]) + np.interp(u, L3, res['o3'][str(g)])
            law[g] = float(np.trapezoid(np.exp(ls) / u * np.exp(u), u)) * norm
        out = {}
        for name, (lo, hi) in [('%d-%d' % b, b) for b in BANDS] + [('pooled', (276, 500))]:
            c = sum(res['counts'].get(str(g), 0) for g in GS if lo <= g <= hi); lw = sum(law[g] for g in GS if lo <= g <= hi)
            out[name] = {'count': c, 'law': lw, 'ratio': c / lw, 'z': (c - lw) / math.sqrt(lw),
                         'z_vs_0.84': (c - 0.84 * lw) / math.sqrt(0.84 * lw), 'g_over_lnx': [lo / math.log(L1), hi / math.log(L0)]}
            print(name, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in out[name].items()}, flush=True)
        res['bands'] = out; res['law_per_gap'] = {str(g): law[g] for g in GS}
        json.dump(res, open(PATH, 'w'), indent=1)
