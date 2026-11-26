# LP for V* = inf F(0) over even F supported in [-1,1] with Fhat(0) = -1 and Fhat(xi) >= 0 for |xi| >= 1.
# Fhat(xi) = int F(x) e^{-2 pi i x xi} dx.  Basis: phi_j(x) = (1-x^2)^p cos(pi j x), j = 0..n.
import numpy as np, sys
from scipy.optimize import linprog
from numpy.polynomial.legendre import leggauss
def run(n, p, Xi=80.0, dxi=0.01):
    xg, wg = leggauss(400)                    # Gauss-Legendre on [-1,1]
    x = 0.5*(xg + 1); w = 0.5*wg              # nodes on [0,1]
    def basis(xx): return np.stack([(1-xx**2)**p*np.cos(np.pi*j*xx) for j in range(n+1)], axis=1)
    B = basis(x)                              # (nodes, n+1)
    def fhat(xi):                             # (len(xi), n+1)
        C = np.cos(2*np.pi*np.outer(xi, x))
        return 2*(C*w) @ B
    xis = np.arange(1.0, Xi, dxi)
    H = fhat(xis); h0 = fhat(np.array([0.0]))[0]
    phi0 = basis(np.array([0.0]))[0]
    r = linprog(phi0, A_ub=-H, b_ub=np.zeros(len(xis)), A_eq=h0[None,:], b_eq=[-1.0],
                bounds=[(None,None)]*(n+1), method='highs')
    if r.status != 0: return None
    c = r.x
    xt = np.arange(Xi, 600, 0.05); tail = (fhat(xt) @ c).min()
    return r.fun, tail
for p in [2, 3]:
    for n in [4, 8, 16, 32, 48]:
        out = run(n, p)
        if out is None: print("p=%d n=%d infeasible/unbounded" % (p, n)); continue
        print("p=%d n=%2d  min F(0) = %.6f   (tail min of Fhat on [80,600] = %.2e)" % (p, n, out[0], out[1])); sys.stdout.flush()
