"""Replacement test of round 31 (non-overlapping intervals in [1e9,3e9)) from my own twin counts."""
import numpy as np, math, json
S2 = 1.3203236316937391
c6 = np.load('tw6_1e9_3e9.npy').astype(float)
law = {10**6: 0.7031, 10**7: 0.6146}       # own law, r31_law.py
out = {}
for H in (10**6, 10**7):
    c = c6 if H == 10**6 else c6.reshape(-1, 10).sum(1)
    m = 10**9 + H * np.arange(len(c)) + H / 2
    e = H * S2 / np.log(m) ** 2
    res = c - e
    K = 10; per = len(c) // K
    meas = [np.var(res[k*per:(k+1)*per], ddof=1) / e[k*per:(k+1)*per].mean() for k in range(K)]
    nodm = [np.mean(res[k*per:(k+1)*per] ** 2) / e[k*per:(k+1)*per].mean() for k in range(K)]
    r = np.array(meas) - law[H]
    out[H] = {'n': len(c), 'mean_resid_over_sd': float(res.mean() / res.std() * math.sqrt(len(c))),
              'meas': float(np.mean(meas)), 'meas_no_demean': float(np.mean(nodm)), 'law': law[H],
              't9': float(r.mean() / (r.std(ddof=1) / math.sqrt(K))), 'se_meas': float(np.std(meas, ddof=1) / math.sqrt(K)),
              'mean_count': float(c.mean()), 'mean_law': float(e.mean())}
    print(H, out[H])
json.dump({str(k): v for k, v in out.items()}, open('r31b_check.json', 'w'), indent=1)
