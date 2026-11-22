"""Round 57 (theory, no data): the entropy of the full Hardy-Littlewood gap law (conditioned expansion, order 2, y = 13)
at fixed lam, against the closed asymptote log2(e lam) - K/ln 2 of round 52, at lam beyond the reach of data.
The distribution is taken over even g up to G = 16 lam (tail mass below 1e-6) and renormalised. usage: r57.py lam1,lam2,..."""
import os, sys, json, math, time
import numpy as np
import r47
from r52 import K

OUT = r47.OUT
PATH = os.path.join(OUT, 'r57.json')
if __name__ == '__main__':
    res = json.load(open(PATH)) if os.path.exists(PATH) else {}
    for lam in map(float, sys.argv[1].split(',')):
        key = '%.1f' % lam if abs(lam - round(lam, 1)) < 1e-9 else '%.4f' % lam
        if key in res:
            continue
        t0 = time.time(); G = int(16 * lam) // 2 * 2
        sh = []
        for g in range(2, G + 1, 2):
            st = r47.setup(g); sh.append(r47.share(st, lam, 2))
        p = np.array(sh); tot = p.sum(); p = p / tot
        H = float(-np.sum(p * np.log2(p)))
        asym = (math.log(math.e * lam) - K) / math.log(2)
        res[key] = {'H_law': H, 'H_asym': asym, 'diff': H - asym, 'lambda_times_diff': lam * (H - asym), 'mass': float(tot), 'G': G,
                    'secs': round(time.time() - t0)}
        json.dump(res, open(PATH, 'w'), indent=1)
        print(key, {k: (round(v, 5) if isinstance(v, float) else v) for k, v in res[key].items()}, flush=True)
