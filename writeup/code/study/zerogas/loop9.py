import sys, numpy as np
def cth(u):
    u = np.asarray(u, float); out = np.empty_like(u)
    a = u < 1; ua = u[a]; out[a] = 1 - (1 + ua)*np.log1p(ua) - np.where(ua < 1, (1 - ua)*np.log1p(-np.minimum(ua, 1 - 1e-15)), 0.0)
    b = ~a; ub = u[b]
    with np.errstate(divide='ignore', invalid='ignore'):
        t = np.where(ub > 1, ub*(ub - 1)*np.log(np.maximum(ub - 1, 1e-300)), 0.0)
    out[b] = 1 + 2*ub*ub*np.log(ub) - t - ub*(ub + 1)*np.log(ub + 1)
    return out
if __name__ == '__main__':
    rng = np.random.default_rng(11); N = int(sys.argv[1]); M = int(sys.argv[2])
    ks = np.arange(1, int(1.5*N) + 1); num = np.zeros((M, len(ks))); den = np.zeros((M, len(ks)))
    for it in range(M):
        Z = (rng.standard_normal((N, N)) + 1j*rng.standard_normal((N, N)))/np.sqrt(2); Q, R = np.linalg.qr(Z)
        x = np.sort(np.angle(np.linalg.eigvals(Q*(np.diag(R)/np.abs(np.diag(R)))[None, :])) % (2*np.pi))*N/(2*np.pi)
        Sz = np.arange(N) - x; Sz = Sz - Sz.mean()
        E = np.exp(-2j*np.pi*np.outer(ks, x)/N); rho = E.sum(axis=1); coef = rho/(2j*np.pi*ks)
        W = 2*np.real(coef[:, None]*np.conj(E))                      # wave of each mode at the eigenvalues, (modes, N)
        num[it] = W @ Sz; den[it] = (W*W).sum(axis=1)
    u = ks/N; edges = np.arange(0, 1.5001, 0.05); grp = np.array_split(np.arange(M), 10); worst = 0
    for i in range(len(edges) - 1):
        sel = (u >= edges[i]) & (u < edges[i + 1])
        c = num[:, sel].sum()/den[:, sel].sum(); e = np.array([num[g][:, sel].sum()/den[g][:, sel].sum() for g in grp]).std(ddof=1)/np.sqrt(10)
        wts = den[:, sel].sum(axis=0); th = float((cth(u[sel])*wts).sum()/wts.sum())
        print('matrices N %d u in [%.2f, %.2f): c = %+.4f +- %.4f | formula %+.4f | z = %+.1f' % (N, edges[i], edges[i + 1], c, e, th, (c - th)/e))
