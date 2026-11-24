import numpy as np, sys, json
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
from readings import *
for test in ('B', 'A'):
    S, full, diag, X, primes, ks = load(test)
    rep = json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/r134_%s.json' % test))
    print('=== %s: reported max|level_full| %.4f, max|level_diag| %.4f (threshold 0.004)' % (test, np.abs(rep['level_full']).max(), np.abs(rep['level_diag']).max()))
    theirS = np.array(rep['sigma2'])
    print('   mean(their Sigma2 - mine) per l:', np.round((theirS - S[:, pidx]).mean(0), 5))
    n = len(S)
    for nm, M in (('full', full), ('diag', diag)):
        R = (S - M)[:, pidx]; m = R.mean(0); C = np.cov(R.T) / n
        c2 = m @ np.linalg.solve(C, m)
        print('   level %s: mean %s ; se %s ; chi2 (16 l, cov) %.1f' % (nm, np.round(m, 4), np.round(np.sqrt(np.diag(C)), 4), c2))
    # rho bootstrap
    rng = np.random.default_rng(1); rr = []
    for b in range(1000):
        i = rng.integers(0, n, n); rr.append(readings(S[i], full[i], diag[i], X[i], primes)['rho'])
    print('   rho %.3f bootstrap se %.3f' % (readings(S, full, diag, X, primes)['rho'], np.std(rr, ddof=1)))
