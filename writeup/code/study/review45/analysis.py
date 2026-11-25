import numpy as np, json, sys
from scipy import stats
SCR = __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/'
P = json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/triple/predictions_kappa3.json'))
R = json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/r138_Fn.json'))
dat = {int(k): v for k, v in json.load(open(SCR + 'data_Fn.json')).items()}
F = P['Fn']; Ds = F['D']; ls = np.array(P['l']); gue = np.array(P['gue_kappa3'])
PR = (2, 3, 5, 7, 11, 13)

def arrays(Dsel):
    idx = [Ds.index(D) for D in Dsel]
    out = {k: np.array([F[k][i] for i in idx]) for k in ('full', 'Q0', 'twopt', 'diagonly', 'base', 'kQ', 'kf3', 'sigma2')}
    out['meanL'] = np.array([F['meanL'][i] for i in idx])
    out['k3'] = np.array([dat[D]['k3'] for D in Dsel]); out['m2'] = np.array([dat[D]['m2'] for D in Dsel])
    out['m1'] = np.array([dat[D]['m1'] for D in Dsel]); out['k3trim'] = np.array([dat[D]['k3trim'] for D in Dsel])
    out['k3l'] = np.array([dat[D]['k3l'] for D in Dsel]); out['nz'] = np.array([dat[D]['nz'] for D in Dsel])
    out['X'] = np.array([[1.] + [float(D % p == 0) for p in PR] for D in Dsel])
    out['D'] = np.array(Dsel)
    return out

def ols(X, Y):
    A = np.linalg.inv(X.T @ X); B = A @ X.T @ Y; res = Y - X @ B
    dof = X.shape[0] - X.shape[1]; Cres = res.T @ res/dof
    return B, A, Cres, res

def readings(a, y=None, alt=None, weights='ols'):
    """a: arrays; y: data kappa3 (n x 16). Returns dict of readings."""
    y = a['k3'] if y is None else y
    full = a['full'] if alt is None else alt
    n = len(y); d = y - full; m = d.mean(0); C = np.cov(d.T)
    T2 = n*m @ np.linalg.solve(C, m); k = 16
    Fst = T2*(n - k)/((n - 1)*k); pT = 1 - stats.f.cdf(Fst, k, n - k)
    aa = full.mean(0) - gue; bb = y.mean(0) - gue
    gA = (aa @ np.linalg.solve(C, bb))/(aa @ np.linalg.solve(C, aa))
    X = a['X']
    B, A, Cres, res = ols(X, y)
    pf = ols(X, full)[0]; pq = ols(X, a['Q0'])[0]; p2 = ols(X, a['twopt'])[0]; pdg = ols(X, a['diagonly'])[0]
    se = np.sqrt(np.outer(np.diag(A), np.diag(Cres)))     # 7 x 16
    rows = [1, 2, 3]
    def amp(palt, w):
        num = sum((w[r]*(B[r] - palt[r])*(pf[r] - palt[r])).sum() for r in rows)
        den = sum((w[r]*(pf[r] - palt[r])**2).sum() for r in rows)
        return num/den
    W = 1/se**2 if weights == 'ols' else np.ones_like(se)
    out = dict(m=m, T2=T2, pT=pT, gA=gA, gQ=amp(pq, W), goff=amp(pdg, W), g2=amp(p2, W), B=B, se=se, pf=pf, pq=pq, pdg=pdg, p2=p2)
    # GLS version of gamma_Q: full covariance of the 48 coefficients: A_rr' Cres
    idx = [(r, j) for r in rows for j in range(16)]
    V = np.array([[A[r, r2]*Cres[j, j2] for (r2, j2) in idx] for (r, j) in idx])
    dQ = np.array([pf[r, j] - pq[r, j] for (r, j) in idx]); obs = np.array([B[r, j] - pq[r, j] for (r, j) in idx])
    Vi = np.linalg.inv(V)
    out['gQ_gls'] = (dQ @ Vi @ obs)/(dQ @ Vi @ dQ); out['gQ_gls_se'] = 1/np.sqrt(dQ @ Vi @ dQ)
    dO = np.array([pf[r, j] - pdg[r, j] for (r, j) in idx]); obsO = np.array([B[r, j] - pdg[r, j] for (r, j) in idx])
    out['goff_gls'] = (dO @ Vi @ obsO)/(dO @ Vi @ dO); out['goff_gls_se'] = 1/np.sqrt(dO @ Vi @ dO)
    # chi2 per prime against law and alternatives
    ch = {}
    for r, p in zip(range(1, 7), PR):
        Vp = A[r, r]*Cres; Vpi = np.linalg.inv(Vp)
        ch[p] = {nm: float((B[r] - pp[r]) @ Vpi @ (B[r] - pp[r])) for nm, pp in (('full', pf), ('Q0', pq), ('diag', pdg), ('GUE', 0*pf))}
    out['chi2'] = ch
    return out

def boot(a, nb=2000, seed=1371, **kw):
    rng = np.random.default_rng(seed); n = len(a['k3']); acc = []
    for b in range(nb):
        ii = rng.integers(0, n, n)
        ab = {k: (v[ii] if isinstance(v, np.ndarray) and v.shape[:1] == (n,) else v) for k, v in a.items()}
        try:
            r = readings(ab, **kw)
            acc.append([r['gA'], r['gQ'], r['goff'], r['g2'], r['gQ_gls'], *r['m']])
        except np.linalg.LinAlgError:
            pass
    acc = np.array(acc)
    return acc.std(0, ddof=1)

if __name__ == '__main__':
    excl = [1153]
    Dsel = [D for D in Ds if D not in excl]
    a = arrays(Dsel)
    # 1. compare own kappa3 with the repo's per-L-function kappa3
    rk = np.array(R['kappa3'])
    print('own kappa3 vs r138_Fn.json per cell: max |diff| = %.2e' % np.max(np.abs(rk - a['k3'])))
    r = readings(a)
    print('m(l) own:', np.round(r['m']*1e4, 2), ' x1e-4'); print('m(l) repo:', np.round(np.array(R['m'])*1e4, 2))
    print('T2 %.3f p %.3f | gA %.4f gQ %.4f goff %.4f g2 %.5f | gQ_gls %.4f +- %.4f  goff_gls %.4f +- %.4f' % (
        r['T2'], r['pT'], r['gA'], r['gQ'], r['goff'], r['g2'], r['gQ_gls'], r['gQ_gls_se'], r['goff_gls'], r['goff_gls_se']))
    print('repo: T2 %.3f p %.3f gA %.4f gQ %.4f goff %.4f g2 %.5f' % (R['T2'], R['T2_p'], R['gamma_A'], R['gamma_Q'], R['gamma_off'], R['gamma_2pt']))
    ru = readings(a, weights='unw'); print('unweighted: gQ %.4f goff %.4f g2 %.5f' % (ru['gQ'], ru['goff'], ru['g2']))
    for p in PR: print('chi2 p=%d' % p, {k: round(v, 1) for k, v in r['chi2'][p].items()})
    se = boot(a, 2000)
    print('bootstrap se: gA %.4f gQ %.4f goff %.4f g2 %.5f gQ_gls %.4f' % tuple(se[:5]))
    print('bootstrap se m(l) x1e-4:', np.round(se[5:]*1e4, 2))
    np.save(SCR + 'boot_se_Fn.npy', se)
