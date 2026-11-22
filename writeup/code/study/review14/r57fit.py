import json, math, numpy as np
d = json.load(open('/home/bas/praime/writeup/data/study/r57.json'))
lam = np.array([float(k) for k in d]); y = np.array([d[k]['diff'] for k in d])
# check asymptote values
K = 0.7627474
for l in lam: print(l, (math.log(math.e*l)-K)/math.log(2), d['%.1f'%l]['H_asym'])
X = np.stack([np.ones_like(lam), 1/lam, np.log(lam)/lam], 1)
c, res, *_ = np.linalg.lstsq(X, y, rcond=None); print('fit diff = c0 + a/lam + b ln(lam)/lam:', c, 'max resid', np.max(np.abs(X@c-y)))
X2 = np.stack([1/lam, np.log(lam)/lam], 1); c2 = np.linalg.lstsq(X2, y, rcond=None)[0]
print('c0=0 fit:', c2, 'max resid', np.max(np.abs(X2@c2-y)), 'zero crossing lam', math.exp(-c2[0]/c2[1]))
for l in (200, 500, 1000, 1e4): print('extrap lam', l, 'diff', (c2[0]+c2[1]*math.log(l))/l)
# three-point fits on subsets (stability of c0)
for idx in ([0,1,2],[3,4,5],[2,3,4,5]):
    cc = np.linalg.lstsq(X[idx], y[idx], rcond=None)[0]; print('subset', lam[idx], 'c0', cc[0])
