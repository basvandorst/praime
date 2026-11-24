import numpy as np, sys
from scipy import stats
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
from readings import *

def coefs(X, Y):
    b, *_ = np.linalg.lstsq(X, Y, rcond=None)
    return b, Y - X @ b

for test in ('B', 'A'):
    S, full, diag, X, primes, ks = load(test)
    Sp = S[:, pidx]; Fp = full[:, pidx]; Dp = diag[:, pidx]
    n, k = X.shape; npr = k - 1
    XtXi = np.linalg.inv(X.T @ X)
    o = readings(S, full, diag, X, primes)
    # OLS se from S alone vs from S - full
    _, eS = coefs(X, Sp); _, eR = coefs(X, Sp - Fp); _, eD = coefs(X, Sp - Dp)
    sS = np.sqrt((eS ** 2).sum(0) / (n - k)); sR = np.sqrt((eR ** 2).sum(0) / (n - k)); sD = np.sqrt((eD ** 2).sum(0) / (n - k))
    print('=== test', test)
    print('residual std of S on X       ', np.round(sS, 4))
    print('residual std of S-full on X  ', np.round(sR, 4))
    print('residual std of S-diag on X  ', np.round(sD, 4))
    # chi2 with proper se (from S - law regression), cells treated independent
    bR, _ = coefs(X, Sp - Fp); bD, _ = coefs(X, Sp - Dp)
    seR = np.sqrt(np.outer(np.diag(XtXi), (eR ** 2).sum(0) / (n - k)))
    seD = np.sqrt(np.outer(np.diag(XtXi), (eD ** 2).sum(0) / (n - k)))
    c2f_proper = ((bR[1:] / seR[1:]) ** 2).sum(); c2d_proper = ((bD[1:] / seD[1:]) ** 2).sum()
    print('registered chi2 (OLS se of S): full %.1f diag %.1f' % (o['chi2f'], o['chi2d']))
    print('chi2 with se of the residual law (diag cov ignored): full %.1f diag %.1f on %d cells' % (c2f_proper, c2d_proper, bR[1:].size))
    # GLS: cov(vec coef) = Sigma_res (16x16) kron XtXi[1:,1:]
    for nm, e, b in (('full', eR, bR), ('diag', eD, bD)):
        Sig = e.T @ e / (n - k)
        V = np.kron(Sig, XtXi[1:, 1:])   # ordering: l-major, prime-minor  -> vec of b[1:].T
        r = b[1:].T.ravel()
        c2 = r @ np.linalg.solve(V, r)
        m = r.size
        # Hotelling-type correction: T2 * (nu - m + 1)/(nu m) ~ F(m, nu - m + 1) ... approximate with nu = n - k
        nu = n - k
        Fst = c2 * (nu - 16 + 1) / (nu * 16) if False else None
        p_chi = stats.chi2.sf(c2, m)
        print('GLS chi2 %s law: %.1f on %d (chi2 p = %.2e)' % (nm, c2, m, p_chi))
    # Null simulation for the registered statistics: synthetic S* = full + resampled residual rows (of S - full on X)
    rng = np.random.default_rng(7)
    sims = []
    for it in range(1000):
        i = rng.integers(0, n, n)
        Sstar = S.copy()
        # residual rows on the full 64 grid
        _, eRall = coefs(X, S - full)
        Sstar = full + eRall[i]
        oo = readings(Sstar, full, diag, X, primes)
        # GLS chi2 of the full law on synthetic data
        bb, ee = coefs(X, Sstar[:, pidx] - Fp)
        Sig = ee.T @ ee / (n - k); V = np.kron(Sig, XtXi[1:, 1:]); r = bb[1:].T.ravel()
        sims.append((oo['chi2f'], oo['gamma'], oo['rho'], r @ np.linalg.solve(V, r)))
    sims = np.array(sims)
    print('NULL (full law true, real residual rows resampled): registered chi2_full mean %.1f, 95%% %.1f, 99%% %.1f; observed %.1f -> p = %.3f' % (
        sims[:, 0].mean(), np.quantile(sims[:, 0], .95), np.quantile(sims[:, 0], .99), o['chi2f'], (sims[:, 0] >= o['chi2f']).mean()))
    print('NULL gamma: mean %.4f sd %.4f; observed %.4f -> two-sided p = %.3f' % (sims[:, 1].mean(), sims[:, 1].std(), o['gamma'], 2 * min((sims[:, 1] <= o['gamma']).mean(), (sims[:, 1] >= o['gamma']).mean())))
    print('NULL rho : mean %.4f sd %.4f; observed %.4f' % (sims[:, 2].mean(), sims[:, 2].std(), o['rho']))
    print('NULL GLS chi2_full: mean %.1f 95%% %.1f 99%% %.1f' % (sims[:, 3].mean(), np.quantile(sims[:, 3], .95), np.quantile(sims[:, 3], .99)))
    # diag law synthetic: S* = diag + residual rows  -> distribution of registered statistics when diag law true
    sims2 = []
    for it in range(300):
        i = rng.integers(0, n, n)
        _, eRall = coefs(X, S - full)
        Sstar = diag + eRall[i]
        oo = readings(Sstar, full, diag, X, primes)
        sims2.append((oo['chi2f'], oo['chi2d'], oo['gamma'], oo['rho']))
    sims2 = np.array(sims2)
    print('IF DIAG LAW TRUE: chi2_full mean %.1f (min %.1f), chi2_diag mean %.1f, gamma mean %.3f sd %.3f, rho mean %.3f' % (
        sims2[:, 0].mean(), sims2[:, 0].min(), sims2[:, 1].mean(), sims2[:, 2].mean(), sims2[:, 2].std(), sims2[:, 3].mean()))
