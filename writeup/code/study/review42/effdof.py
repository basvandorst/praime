import numpy as np, sys
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
from readings import *
for test in ('B','A'):
    S, full, diag, X, primes, ks = load(test)
    Y = (S - full)[:, pidx]
    b, *_ = np.linalg.lstsq(X, Y, rcond=None); e = Y - X @ b
    C = np.corrcoef(e.T)
    lam = np.linalg.eigvalsh(C)
    XtXi = np.linalg.inv(X.T @ X); D = np.sqrt(np.diag(XtXi)); Cp = XtXi / np.outer(D, D)
    lp = np.linalg.eigvalsh(Cp[1:, 1:])
    print(test, 'corr(l=6,6.5)=%.2f corr(l=7,8)=%.2f corr(l=0.5,8)=%.2f' % (C[11, 12], C[13, 15], C[0, 15]),
          '| effective independent lengths (sum lam)^2/sum lam^2 = %.1f of 16' % (lam.sum() ** 2 / (lam ** 2).sum()),
          '| max |corr| between prime coefficients %.2f' % np.abs(Cp[1:, 1:] - np.eye(len(primes))).max())
