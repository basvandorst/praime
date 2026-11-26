# V* = inf F(0) = inf sum_k G(k), G = Fhat in PW (type 2 pi, F supported in [-1,1], F even),
# G(0) = -1, G >= 0 on |xi| >= 1.  Hermite basis: G(xi) = sum_{|k|<=K} [g_k + (xi-k) d_k] sinc^2(pi(xi-k)),
# g even, d odd.  Poisson: F(0) = sum_k g_k.
import numpy as np, sys
from scipy.optimize import linprog
def sinc2(z): return np.sinc(z)**2           # np.sinc(z) = sin(pi z)/(pi z)
def run(K, dxi=0.005, far=40):
    ks = np.arange(1, K+1)
    xis = np.concatenate([np.arange(1.0, K + far, dxi), np.arange(K+far, 20*K+400, 0.05)])
    # columns: g_0, g_1..g_K, d_1..d_K
    cols = [sinc2(xis)]
    for k in ks: cols.append(sinc2(xis-k) + sinc2(xis+k))
    for k in ks: cols.append((xis-k)*sinc2(xis-k) - (xis+k)*sinc2(xis+k))
    M = np.stack(cols, axis=1)
    # asymptotic constraint: sum_k g_k + 2 sum k d_k >= 0  (coefficient of sin^2/(pi xi)^2)
    asym = np.concatenate([[1.0], 2*np.ones(K), 2*ks.astype(float)])
    A_ub = -np.vstack([M, asym[None,:]]); b_ub = np.zeros(A_ub.shape[0])
    c = np.concatenate([[1.0], 2*np.ones(K), np.zeros(K)])   # F(0) = g_0 + 2 sum g_k
    A_eq = np.zeros((1, 2*K+1)); A_eq[0,0] = 1.0
    bounds = [(None,None)] + [(0,None)]*K + [(None,None)]*K
    r = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=[-1.0], bounds=bounds, method='highs')
    return r
for K in [2, 4, 8, 16, 32, 64]:
    r = run(K)
    if r.status != 0: print("K=%d status %d" % (K, r.status)); continue
    g = r.x[1:K+1]; d = r.x[K+1:]
    print("K=%2d  V*(K) = min F(0) = %.6f   sum_{k>=1} g_k = %.5f  g_1..3 = %s" % (K, r.fun, g.sum(), np.round(g[:3],5))); sys.stdout.flush()
