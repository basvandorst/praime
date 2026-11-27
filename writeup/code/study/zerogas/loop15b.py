import numpy as np
from loop15 import c_exact
rng = np.random.default_rng(1515); N, M = 10, 200000
ks = np.arange(1, N + 1); num = np.zeros((M, N)); den = np.zeros((M, N))
for it in range(M):
    Z = (rng.standard_normal((N, N)) + 1j*rng.standard_normal((N, N)))/np.sqrt(2); Q, R = np.linalg.qr(Z)
    x = np.sort(np.angle(np.linalg.eigvals(Q*(np.diag(R)/np.abs(np.diag(R)))[None, :])) % (2*np.pi))*N/(2*np.pi)
    Sz = np.arange(N) - x
    E = np.exp(-2j*np.pi*np.outer(ks, x)/N); coef = E.sum(axis=1)/(2j*np.pi*ks); W = 2*np.real(coef[:, None]*np.conj(E))
    num[it] = W @ Sz; den[it] = (W*W).sum(axis=1)
grp = np.array_split(np.arange(M), 40)
for i, k in enumerate(ks):
    c = num[:, i].mean()/den[:, i].mean(); e = np.array([num[g, i].mean()/den[g, i].mean() for g in grp]).std(ddof=1)/np.sqrt(40)
    ex = float(c_exact(N, int(k)))
    # the two expectations separately, against their exact values (units 1/(2 pi^2 k))
    print('N %d k %2d: simulation %+.5f +- %.5f | exact %+.5f | z = %+.1f | E num, E den in units N/(2 pi^2 k): %.5f, %.5f' % (N, k, c, e, ex, (c - ex)/e, num[:, i].mean()*2*np.pi**2*k/N, den[:, i].mean()*2*np.pi**2*k/N))
