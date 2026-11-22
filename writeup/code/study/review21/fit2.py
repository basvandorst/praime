import json, numpy as np
d = json.load(open('/home/bas/praime/writeup/data/study/r70.json'))
L = [float(k) for k in d if k != 'fit']; Y = [d[k]['lambda_times_diff'] for k in d if k != 'fit']
L = np.array(L+[600.]); Y = np.array(Y+[1.3759788846486543]); l = np.log(L)
for name, X in [('a+b ln', np.c_[np.ones_like(l), l]), ('a+b ln+c ln2', np.c_[np.ones_like(l), l, l*l]), ('a+b ln + c/lam', np.c_[np.ones_like(l), l, 1/L])]:
    c, *_ = np.linalg.lstsq(X, Y, rcond=None); print(name, c.round(4), 'maxres', abs(Y-X@c).max().round(4))
print('local slopes', (np.diff(Y)/np.diff(l)).round(3))
