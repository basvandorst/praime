"""Round 30 (preregister/round30.txt): counts of large gaps in [1e9, 2e9) against the law + overlap correction."""
import os, json, math, pickle
import numpy as np
import theory2 as T2
import theory_r3 as R3
import common as C

OUT = R3.OUT
cache = pickle.load(open(R3.CACHE, 'rb'))
L, W = 10 ** 9, 10 ** 9
BANDS = [(120, 178), (180, 238), (240, 298), (300, 358), (360, 454)]
hist = np.zeros(2000, dtype=np.int64); nprimes = 0; last = None
for s in range(100):
    lo = L + s * 10 ** 7
    P = np.nonzero(C.primes_mask(lo, 10 ** 7))[0].astype(np.int64) + lo
    if last is not None:
        P = np.concatenate(([last], P))
    g = np.diff(P)
    hist += np.bincount(np.minimum(g, 1999), minlength=2000)
    nprimes += len(P) - (1 if last is not None else 0); last = P[-1]
ngaps = int(hist.sum())
u = np.linspace(math.log(L), math.log(L + W), 201); w = np.exp(u) / u
share = {}
for h in range(2, R3.HMAX + 1, 2):
    SB, la = cache[(0, h)]
    if SB == 0:
        share[h] = 0.0; continue
    pt = T2.Pattern([0, h], 2)
    corr = np.array([T2.overlap_correction(pt, l) for l in R3.LAMG])
    share[h] = float(np.trapezoid(w * SB / u * np.exp(np.interp(u, R3.LAMG, la + corr)), u) / np.trapezoid(w, u))
tot = sum(share.values())
out = {'gaps': ngaps, 'law_norm': tot, 'bands': {}}
lnx = math.log(1.5e9)
for a, b in BANDS:
    cnt = int(hist[a:b + 1:2].sum()) if a % 2 == 0 else int(hist[a:b + 1].sum())
    cnt = int(sum(hist[h] for h in range(a, b + 1, 2)))
    law = sum(share[h] for h in range(a, b + 1, 2)) / tot * ngaps
    out['bands']['%d-%d' % (a, b)] = {'g_over_lnx': [a / lnx, b / lnx], 'count': cnt, 'law': law, 'ratio': cnt / law,
                                      'z': (cnt - law) / math.sqrt(law)}
    print(a, b, '%.1f-%.1f' % (a / lnx, b / lnx), cnt, '%.1f' % law, 'ratio %.3f' % (cnt / law), 'z %.1f' % ((cnt - law) / math.sqrt(law)), flush=True)
json.dump(out, open(os.path.join(OUT, 'r30.json'), 'w'), indent=1)
