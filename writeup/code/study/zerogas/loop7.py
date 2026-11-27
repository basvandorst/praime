import sys, numpy as np
rng = np.random.default_rng(7)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 400; M = int(sys.argv[2]) if len(sys.argv) > 2 else 300
edges = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2]
nb = len(edges) - 1
num_c = np.zeros((M, nb)); den_c = np.zeros((M, nb)); num_b = np.zeros((M, nb)); den_b = np.zeros((M, nb))
ks = np.arange(1, int(1.2*N) + 1); band = np.digitize(ks/N, edges) - 1
for it in range(M):
    Z = (rng.standard_normal((N, N)) + 1j*rng.standard_normal((N, N)))/np.sqrt(2); Q, R = np.linalg.qr(Z)
    U = Q*(np.diag(R)/np.abs(np.diag(R)))[None, :]
    x = np.sort(np.angle(np.linalg.eigvals(U)) % (2*np.pi))*N/(2*np.pi)
    Sz = np.arange(N) - x; Sz = Sz - Sz.mean()                                  # S at the eigenvalues (mean of both sides), up to a constant
    xm = (x + np.roll(x, -1))/2; xm[-1] += N/2; Sm = (Sz + np.roll(Sz, -1))/2   # midpoints (circular)
    rho = np.exp(-2j*np.pi*np.outer(ks, x)/N).sum(axis=1)                       # modes of the density
    coef = rho/(2j*np.pi*ks)
    for b in range(nb):
        sel = band == b
        if not sel.any(): continue
        Wz = 2*np.real((coef[sel, None]*np.exp(2j*np.pi*np.outer(ks[sel], x)/N)).sum(axis=0))
        Wm = 2*np.real((coef[sel, None]*np.exp(2j*np.pi*np.outer(ks[sel], xm)/N)).sum(axis=0))
        num_c[it, b] = Sz @ Wz; den_c[it, b] = Wz @ Wz; num_b[it, b] = Sm @ Wm; den_b[it, b] = Wm @ Wm
def stat(num, den):
    v = num.sum(0)/den.sum(0); g = np.array_split(np.arange(M), 10)
    e = np.array([num[i].sum(0)/den[i].sum(0) for i in g]).std(axis=0, ddof=1)/np.sqrt(10); return v, e
c, ce = stat(num_c, den_c); b, be = stat(num_b, den_b)
for i in range(nb):
    print('CUE N %d u in [%.2f, %.2f): c(u) = %+.4f +- %.4f   b(u) = %+.4f +- %.4f' % (N, edges[i], edges[i + 1], c[i], ce[i], b[i], be[i]))
np.save('cue_N%d.npy' % N, np.array([c, ce, b, be]))
