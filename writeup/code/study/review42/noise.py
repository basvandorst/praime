import numpy as np, sys
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
from readings import *
for test in ('B', 'A'):
    S, full, diag, X, primes, ks = load(test)
    Sp = S[:, pidx]; Fp = full[:, pidx]
    r = Sp - Fp
    print(test, 'std over L-functions of (S - full) at l=0.5..8:', np.round(r.std(0, ddof=1), 4))
    print(test, 'std of full prediction itself      :', np.round(Fp.std(0, ddof=1), 4))
