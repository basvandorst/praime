import json, numpy as np
d = json.load(open('/home/bas/praime/writeup/data/study/r70.json'))
L = np.array([float(k) for k in d if k != 'fit']); Y = np.array([d[k]['lambda_times_diff'] for k in d if k != 'fit'])
X = np.c_[np.ones_like(L), np.log(L)]
c, *_ = np.linalg.lstsq(X, Y, rcond=None); r = Y - X@c
print('fit a,b', c, 'maxres', abs(r).max(), r.round(4))
print('pred 600', c[0]+c[1]*np.log(600), 'computed 1.37598')
X3 = np.c_[np.ones_like(L), np.log(L), 1/L]
c3, *_ = np.linalg.lstsq(X3, Y, rcond=None); print('a+b ln + c/lam', c3, abs(Y-X3@c3).max())
for lo in [21, 60]:
    m = L >= lo; c2, *_ = np.linalg.lstsq(X[m], Y[m], rcond=None); print('lam>=', lo, c2)
# pure a/lam on data-range lambdas
m = L <= 30; print('mean lam*diff in data range', Y[m])
