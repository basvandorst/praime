import numpy as np, sys
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
from readings import *
S, full, diag, X, primes, ks = load('B')
n = len(S)
b, *_ = np.linalg.lstsq(X, S - full, rcond=None); E = (S - full) - X @ b
rng = np.random.default_rng(11)
for a in (0.6, 0.7, 0.75, 0.8, 0.9, 1.0, 1.1, 1.2, 1.25, 1.3, 1.4):
    ok = 0; N = 200
    for it in range(N):
        Ss = diag + a * (full - diag) + E[rng.integers(0, n, n)]
        o = readings(Ss, full, diag, X, primes)
        p4 = np.all(np.sign(o['p2']) == np.array([-1, 1, -1, 1, -1]))
        if 0.7 <= o['gamma'] <= 1.3 and o['chi2f'] <= 96 and 0.40 <= o['rho'] <= 0.72 and p4: ok += 1
    print('true amplitude %.2f x BK off: all four confirmation readings pass in %3.0f%% of simulated data sets' % (a, 100 * ok / N))
