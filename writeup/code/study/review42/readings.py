import numpy as np, json, sys
W = __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/'
d = json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/bk/predictions_per_lfunction.json'))
lgrid = np.array(d['l'])
prim = np.arange(0.5, 8.01, 0.5)
pidx = np.array([int(np.argmin(abs(lgrid - l))) for l in prim])

def load(test):
    if test == 'B':
        z = np.load(W + 'S_B.npz'); S = z['S']; ks = [abs(x) for x in z['D']]; primes = [2, 3, 5, 7, 11, 13]
        full = np.array(d['B']['full']); diag = np.array(d['B']['diag'])
    else:
        z = np.load(W + 'S_A.npz'); S = z['S']; ks = list(z['N']); primes = [2, 3, 5, 7]
        full = np.array(d['A']['full']); diag = np.array(d['A']['diag'])
    X = np.column_stack([np.ones(len(ks))] + [[1.0 if k % p == 0 else 0.0 for k in ks] for p in primes])
    return S, full, diag, X, primes, ks

def ols(X, Y):
    # Y: n x m ; returns coef (p x m), se (p x m)
    beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
    resid = Y - X @ beta
    n, p = X.shape
    s2 = (resid ** 2).sum(0) / (n - p)
    XtXi = np.linalg.inv(X.T @ X)
    se = np.sqrt(np.outer(np.diag(XtXi), s2))
    return beta, se

def readings(S, full, diag, X, primes, verbose=True):
    Sp = S[:, pidx]; Fp = full[:, pidx]; Dp = diag[:, pidx]
    bf, se = ols(X, Sp); bF, _ = ols(X, Fp); bD, _ = ols(X, Dp)
    fit = bf[1:]; se = se[1:]; pf = bF[1:]; pd = bD[1:]
    off = pf - pd; res = fit - pd
    gamma = (res * off).sum() / (off ** 2).sum()
    chi2f = (((fit - pf) / se) ** 2).sum(); chi2d = (((fit - pd) / se) ** 2).sum()
    rho = (fit[:, 0] * pd[:, 0]).sum() / (pd[:, 0] ** 2).sum()
    rho_pred = (pf[:, 0] * pd[:, 0]).sum() / (pd[:, 0] ** 2).sum()
    p2 = res[0, 11:16]; p2p = off[0, 11:16]; p2se = se[0, 11:16]
    lev = (Sp - Fp).mean(0); levd = (Sp - Dp).mean(0)
    # pooled per-L-function: (S - diag) on (full - diag) with free intercept per l
    y = (Sp - Dp); x = (Fp - Dp)
    yc = y - y.mean(0); xc = x - x.mean(0)
    pooled = (yc * xc).sum() / (xc ** 2).sum()
    out = dict(gamma=gamma, chi2f=chi2f, chi2d=chi2d, rho=rho, rho_pred=rho_pred, p2=p2, p2p=p2p, p2se=p2se,
               lev=lev, levd=levd, pooled=pooled, fit=fit, se=se, pf=pf, pd=pd)
    return out

def boot(S, full, diag, X, primes, nb=2000, seed=134):
    rng = np.random.default_rng(seed)
    n = len(S); g = []; pl = []
    for b in range(nb):
        i = rng.integers(0, n, n)
        try:
            o = readings(S[i], full[i], diag[i], X[i], primes)
        except np.linalg.LinAlgError:
            continue
        g.append(o['gamma']); pl.append(o['pooled'])
    g = np.array(g); pl = np.array(pl)
    return g.std(ddof=1), np.quantile(g, [0.025, 0.975]), pl.std(ddof=1)

if __name__ == '__main__':
    test = sys.argv[1]
    S, full, diag, X, primes, ks = load(test)
    # compare with their sigma2
    rep = json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/r134_%s.json' % test))
    theirS = np.array(rep['sigma2'])
    print('their sigma2 shape', theirS.shape)
    if theirS.shape[1] == 16:
        print('max |mine - theirs| (primary grid)', np.abs(S[:, pidx] - theirS).max())
    else:
        print('max |mine - theirs|', np.abs(S - theirS).max())
    o = readings(S, full, diag, X, primes)
    print('gamma %.4f chi2 full %.1f diag %.1f rho %.4f (pred %.4f) pooled %.4f' % (o['gamma'], o['chi2f'], o['chi2d'], o['rho'], o['rho_pred'], o['pooled']))
    print('p2 res', np.round(o['p2'], 4), 'pred', np.round(o['p2p'], 4), 'se', np.round(o['p2se'], 4))
    print('level full', np.round(o['lev'], 4))
    print('max |level full|', np.abs(o['lev']).max(), ' max |level diag|', np.abs(o['levd']).max())
    gs, q, ps = boot(S, full, diag, X, primes)
    print('boot gamma se %.4f q %.4f %.4f pooled se %.4f' % (gs, q[0], q[1], ps))
