"""Round 53, final numerics (after the twelfth review): order 2 evaluated directly at Gauss-Legendre nodes in u = ln t
(no interpolation), the order-3 + overlap correction at lam = 26, 30, 34, 38, 42.84 with cubic interpolation.
usage: r53g.py o2 | o3 | compare   (resumable)"""
import os, sys, json, math, time
import numpy as np
from scipy.interpolate import CubicSpline
import common as C
import r47
from r51 import o3_delta

OUT = r47.OUT
PATH = os.path.join(OUT, 'r53g.json')
GS = [40, 50, 60, 70, 72, 100, 150, 200, 250, 300, 350, 400, 454, 460, 482, 500, 520, 538, 556, 580, 610, 640, 670, 700]
A, B = math.log(1e7), math.log(4e18)
XN, WN = np.polynomial.legendre.leggauss(160)
UN = 0.5 * (B - A) * XN + 0.5 * (B + A); WU = 0.5 * (B - A) * WN
L3 = np.array([26.0, 30.0, 34.0, 38.0, B])

if __name__ == '__main__':
    res = json.load(open(PATH)) if os.path.exists(PATH) else {'o2': {}, 'o3': {}}
    mode = sys.argv[1]; t0 = time.time()
    if mode in ('o2', 'o3'):
        for g in GS:
            if str(g) in res[mode] or time.time() - t0 > 420:
                continue
            if mode == 'o2':
                st = r47.setup(g); res['o2'][str(g)] = [math.log(r47.share(st, float(u), 2)) for u in UN]
            else:
                res['o3'][str(g)] = o3_delta(g, L3)
            json.dump(res, open(PATH, 'w'))
            print(mode, g, round(time.time() - t0), flush=True)
        print(mode, 'done', sum(str(g) in res[mode] for g in GS), 'of', len(GS), flush=True)
    else:
        data = json.load(open(os.path.join(OUT, 'oes_gap_counts_4e18.json')))['N']
        M = C.primes_mask(2, 10 ** 7 - 2); P = np.append(np.nonzero(M)[0] + 2, 10000019); g0 = np.diff(P)
        out = {}
        for g in GS:
            base = int(np.sum(g0 == g)); d = data[str(g)]
            ls2 = np.array(res['o2'][str(g)])
            d3 = CubicSpline(L3, res['o3'][str(g)])(np.clip(UN, L3[0], L3[-1]))
            # below lam = 26 the correction is held at its lam = 26 value (weight there is below 1e-7)
            law2 = float(np.sum(WU * np.exp(ls2) * np.exp(UN) / UN)) + base
            law3 = float(np.sum(WU * np.exp(ls2 + d3) * np.exp(UN) / UN)) + base
            out[str(g)] = {'data': d, 'rel_o2': (d - law2) / d, 'rel_o3': (d - law3) / d, 'step': (law3 - law2) / d,
                           'sd_rel': 1 / math.sqrt(d), 'remainder_over_step': abs((d - law3) / (law3 - law2))}
            print(g, ' '.join('%s %+.4f%%' % (k, 100 * out[str(g)][k]) for k in ('rel_o2', 'rel_o3', 'step')), 'rem/step %.2f' % out[str(g)]['remainder_over_step'], flush=True)
        res['compare'] = out
        json.dump(res, open(PATH, 'w'), indent=1)
