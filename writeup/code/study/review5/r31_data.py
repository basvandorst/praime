"""Own recount of the sliding twin-count Var/mean at train for H = 1e4, 1e5 (and 1e3), detrended by the law's mean,
10 blocks of 9e6, as described in STUDY.md round 31; plus non-overlapping H=1e5 intervals as a cross-check."""
import math, json
import numpy as np
from lib import seg_mask

S2 = 1.3203236316937391
L, W = 10 ** 7, 9 * 10 ** 7
HS = (1000, 10000, 100000)
law = {1000: 0.8556, 10000: 0.7752, 100000: 0.6782}   # from r31_law.py (own)
res = {H: [] for H in HS}; nonover = []
for b in range(10):
    lo, n = L + b * W // 10, W // 10
    M = seg_mask(lo, n + max(HS) + 2)
    tw = (M[:-2] & M[2:]).astype(np.int32)
    del M
    c = np.concatenate(([0], np.cumsum(tw)))
    for H in HS:
        cnt = (c[H:n + 1] - c[:n + 1 - H]).astype(np.float64)
        pos = lo + np.arange(len(cnt)) + H / 2
        e = H * S2 / np.log(pos) ** 2
        res[H].append(float(np.var(cnt - e) / e.mean()))
    # non-overlapping H=1e5 intervals, residuals about the law mean, no demeaning
    k = n // 100000
    blk = c[100000 * np.arange(1, k + 1)] - c[100000 * np.arange(k)]
    mid = lo + 100000 * np.arange(k) + 5e4
    e = 1e5 * S2 / np.log(mid) ** 2
    nonover.append((blk - e).tolist())
    print(b, [round(res[H][-1], 4) for H in HS], flush=True)
out = {}
for H in HS:
    r = np.array(res[H]) - law[H]
    t9 = r.mean() / (r.std(ddof=1) / math.sqrt(10))
    out[H] = {'meas': float(np.mean(res[H])), 'law': law[H], 't9_vs_const_law': float(t9)}
    print(H, out[H])
z = np.concatenate([np.array(x) for x in nonover])
e = 1e5 * S2 / math.log(5.5e7) ** 2
out['nonoverlap_1e5'] = {'n': len(z), 'v_about_law_mean': float(np.mean(z ** 2) / e), 'v_demeaned': float(np.var(z, ddof=1) / e)}
print(out['nonoverlap_1e5'])
json.dump(out, open('/home/bas/scratch/study-review5/r31_data.json', 'w'), indent=1)
