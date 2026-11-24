import numpy as np, sys
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
from readings import *

def gam_parts(S, full, diag, X):
    Sp = S[:, pidx]; Fp = full[:, pidx]; Dp = diag[:, pidx]
    bf, _ = ols(X, Sp); bF, _ = ols(X, Fp); bD, _ = ols(X, Dp)
    off = (bF - bD)[1:]; res = (bf - bD)[1:]
    return res, off

for test in ('B', 'A'):
    S, full, diag, X, primes, ks = load(test)
    n = len(S)
    res, off = gam_parts(S, full, diag, X)
    g = (res * off).sum() / (off ** 2).sum()
    rng = np.random.default_rng(134)
    B = []
    for b in range(1000):
        i = rng.integers(0, n, n)
        r_, o_ = gam_parts(S[i], full[i], diag[i], X[i])
        B.append((r_, o_))
    def bs(fun):
        v = np.array([fun(r_, o_) for r_, o_ in B]); return v.std(ddof=1)
    print('=== test', test, 'gamma %.4f' % g)
    print('share of sum off^2 per prime:', {p: round(float((off[j] ** 2).sum() / (off ** 2).sum()), 3) for j, p in enumerate(primes)})
    for j, p in enumerate(primes):
        gp = (res[j] * off[j]).sum() / (off[j] ** 2).sum()
        se = bs(lambda r_, o_, j=j: (r_[j] * o_[j]).sum() / (o_[j] ** 2).sum())
        print('  p=%2d gamma_p %.3f +- %.3f' % (p, gp, se))
    li = np.arange(16)
    for nm, m in (('half-integer l', (prim % 1) == 0.5), ('integer l', (prim % 1) == 0), ('l<=4', prim <= 4), ('l>4', prim > 4), ('l<=2', prim <= 2)):
        gp = (res[:, m] * off[:, m]).sum() / (off[:, m] ** 2).sum()
        se = bs(lambda r_, o_, m=m: (r_[:, m] * o_[:, m]).sum() / (o_[:, m] ** 2).sum())
        print('  %-15s gamma %.3f +- %.3f' % (nm, gp, se))
    gx = (res[1:] * off[1:]).sum() / (off[1:] ** 2).sum()
    print('  without p=2: gamma %.3f +- %.3f' % (gx, bs(lambda r_, o_: (r_[1:] * o_[1:]).sum() / (o_[1:] ** 2).sum())))
    # residual relative to full law per cell, in se units (bootstrap se of res-off)
    d = res - off
    sd = np.array([r_ - o_ for r_, o_ in B]).std(0, ddof=1)
    z = d / sd
    print('  z = (fit - full)/boot se, by prime (rows) and l (cols 0.5..8):')
    for j, p in enumerate(primes):
        print('   p=%2d ' % p + ' '.join('%+5.1f' % v for v in z[j]))
    print('  sum z^2 = %.1f on %d cells' % ((z ** 2).sum(), z.size))
