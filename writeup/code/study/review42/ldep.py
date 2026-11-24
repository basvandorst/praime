import numpy as np, sys
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
from readings import *
S, full, diag, X, primes, ks = load('B')
Ds = np.load(W + 'S_B.npz')['D']
Y = (S - full)[:, pidx]
for nm, v in (('ln|D|', np.log(np.abs(Ds))), ('D<0', (Ds < 0).astype(float)), ('#bad primes', np.array([len(__import__('mypred').badprimes(int(D))) for D in Ds], float))):
    X2 = np.column_stack([X, v - v.mean()])
    b, se = ols(X2, Y)
    print('%-12s coef per l: %s' % (nm, ' '.join('%+.1f' % z for z in b[-1] / se[-1])), ' (t-values)')
