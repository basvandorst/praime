import numpy as np, sys
sys.path.insert(0, '.')
from myfit import load, design, ols
from stress import cp_pred, cl_pred
rows = [r for r in load() if r['D'] not in (197, -395)]
X, nn = design(rows)
f = lambda a: ' '.join('%+.5f' % x for x in ols(X, np.array(a))[nn:nn+4])
print('frozen 3pt  ', f([r['pred']['cum2_3pt'] for r in rows]))
print('compact add.', f([cp_pred(r) for r in rows]))
print('closure     ', f([cl_pred(r) for r in rows]))
nmulti = sum(1 for r in rows if sum(r['q'] % p == 0 for p in (2,3,5,7,11,13)) >= 2)
print('L-functions with >= 2 of the six primes:', nmulti, 'of', len(rows))
