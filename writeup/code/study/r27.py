"""Round 27 (preregister/round27.txt): tails of gaps >= 6, 8, 10 ln x in the windows against the law + overlap correction."""
import os, json, math, pickle
import numpy as np
import theory2 as T2
import theory_r3 as R3
import common as C

OUT = R3.OUT
cache = pickle.load(open(R3.CACHE, 'rb'))
corr = {}


def share(h, L, W):
    SB, la = cache[(0, h)]
    if SB == 0:
        return 0.0
    if h not in corr:
        pt = T2.Pattern([0, h], 2)
        corr[h] = np.array([T2.overlap_correction(pt, l) for l in R3.LAMG])
    u = np.linspace(math.log(L), math.log(L + W), 401); w = np.exp(u) / u
    f = SB / u * np.exp(np.interp(u, R3.LAMG, la + corr[h]))
    return float(np.trapezoid(w * f, u) / np.trapezoid(w, u))


out = {}
for win, (L, W) in R3.WIN.items():
    lm = math.log(L + W / 2)
    ths = [int(round(k * lm / 2)) * 2 for k in (6, 8, 10)]
    gl = {h: share(h, L, W) for h in range(2, R3.HMAX + 1, 2)}
    tot = sum(gl.values())
    M = C.primes_mask(L, W); g = np.diff(np.nonzero(M)[0])
    n = len(g)
    for k, th in zip((6, 8, 10), ths):
        law = sum(v for h, v in gl.items() if h >= th) / tot
        cnt = int(np.sum(g >= th))
        out['%s_%dlnx' % (win, k)] = {'threshold': th, 'count': cnt, 'expected': law * n, 'ratio': cnt / (law * n),
                                      'z_binom': (cnt - law * n) / math.sqrt(law * n)}
        print(win, k, th, cnt, '%.1f' % (law * n), 'ratio %.3f' % (cnt / (law * n)), 'z %.1f' % out['%s_%dlnx' % (win, k)]['z_binom'], flush=True)
json.dump(out, open(os.path.join(OUT, 'r27.json'), 'w'), indent=1)
