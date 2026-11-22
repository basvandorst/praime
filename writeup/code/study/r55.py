"""Round 55 (preregister/round55.txt): gaps of 208..300 (g/ln x >= 10) in [1e8, 1e9) against the conditioned law.
usage: r55.py"""
import os, json, math
import numpy as np
import mpmath
import common as C
import r47
from r51 import o3_delta

L0, L1 = 10 ** 8, 10 ** 9
A, B = math.log(L0), math.log(L1)
XN, WN = np.polynomial.legendre.leggauss(32)
UN = 0.5 * (B - A) * XN + 0.5 * (B + A); WU = 0.5 * (B - A) * WN
L3 = np.array([A, 0.5 * (A + B), B])
GS = list(range(208, 301, 2))
BANDS = [(208, 238), (240, 268), (270, 300)]

if __name__ == '__main__':
    cnt = {}; last = None; ng = 0
    for s in range(L0, L1, 10 ** 8):
        P = np.nonzero(C.primes_mask(s, 10 ** 8))[0].astype(np.int64) + s
        if last is not None:
            P = np.concatenate([[last], P])
        last = int(P[-1]); g = np.diff(P); ng += len(g)
        v, c = np.unique(g[g >= 200], return_counts=True)
        for a, b in zip(v.tolist(), c.tolist()):
            cnt[a] = cnt.get(a, 0) + b
    norm = ng / float(mpmath.li(L1) - mpmath.li(L0))
    law = {}
    for g in GS:
        st = r47.setup(g)
        ls2 = np.array([math.log(r47.share(st, float(u), 2)) for u in UN])
        d3 = np.polyval(np.polyfit(L3, o3_delta(g, L3), 2), UN)
        law[g] = float(np.sum(WU * np.exp(ls2 + d3) * np.exp(UN) / UN)) * norm
    out = {'gaps': ng, 'counts': {str(k): v for k, v in cnt.items()}, 'law': {str(g): law[g] for g in GS}, 'bands': {}}
    for name, (lo, hi) in [('%d-%d' % b, b) for b in BANDS] + [('pooled', (208, 300))]:
        c = sum(cnt.get(g, 0) for g in GS if lo <= g <= hi); lw = sum(law[g] for g in GS if lo <= g <= hi)
        out['bands'][name] = {'count': c, 'law': lw, 'ratio': c / lw, 'z': (c - lw) / math.sqrt(lw),
                              'z_vs_0.54': (c - 0.54 * lw) / math.sqrt(0.54 * lw), 'g_over_lnx': [lo / B, hi / A]}
        print(name, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in out['bands'][name].items()}, flush=True)
    json.dump(out, open(os.path.join(r47.OUT, 'r55.json'), 'w'), indent=1)
