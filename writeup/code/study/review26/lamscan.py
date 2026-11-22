"""Review 26: single-lam law against the 1e7-1e8 data: which lam pass (max|z| < 2, < 4)?"""
import json, math
import numpy as np
from scipy.interpolate import CubicSpline
C = json.load(open('/home/bas/scratch/study-review26/count.json'))
Z = np.load('/home/bas/scratch/study-review26/law_win.npz'); dist = Z['dist']; lams = Z['lams']; Ds = Z['D']
cs = [CubicSpline(lams, dist[:, k, :], axis=0) for k in range(4)]
res = []
for lam in np.linspace(lams[0], lams[-1], 461):
    z = []
    for k in range(4):
        d = cs[k](lam); n = C['data'][str(k + 1)]['n']
        for q in (3, 5, 7, 11, 13):
            m = C['data'][str(k + 1)][str(q)] / n; z.append((d[Ds % q == 0].sum() / d.sum() - m) / math.sqrt(m * (1 - m) / n))
    res.append((lam, np.abs(z).max()))
r = np.array(res)
for t in (2, 4):
    ok = r[r[:, 1] < t, 0]; print('max|z| < %d for lam in [%.3f, %.3f]' % (t, ok.min(), ok.max()))
print('best single lam %.3f max|z| %.2f' % tuple(r[r[:, 1].argmin()]))
