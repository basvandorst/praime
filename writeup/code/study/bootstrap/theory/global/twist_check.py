# Test of the twisted summation formula on a non-self-dual example:
#   zeta_P(s) = zeta(s) * Phi(s),  Phi(s) = 1 + eps/(s-rho) + conj(eps)/(s-conj(rho)),  rho = first zeta zero.
#   Lambda_P(1-s) = B(s) Lambda_P(s) with B = Phi(1-s)/Phi(s) (bounded on Re s >= 1.1 for small eps).
# Mellin measure on [1,inf): sum_n delta_n + D(x) dx,  D(x) = 2 Re( eps x^{rho-1} sum_{n<x} n^{-rho} ).
# Claim:  sum_{n>=1} Fhat(n) + int_1^inf Fhat(x) D(x) dx = (r/2)(F(0) - B(1) Fhat(0)),  r = Phi(1), B(1) = Phi(0)/Phi(1),
# for even F supported in (-1,1).
import numpy as np, mpmath as mp
from numpy.polynomial.legendre import leggauss
rho = complex(0.5, float(mp.im(mp.zetazero(1))))
def Phi(s, eps): return 1 + eps/(s-rho) + np.conj(eps)/(s-np.conj(rho))
p = 8
xg, wg = leggauss(300); xq = 0.5*(xg+1); wq = 0.5*wg
Fq = (1 - xq**2)**p
from scipy.special import jv, gamma
def Fhat(xi):                       # closed form: sqrt(pi) Gamma(p+1) J_{p+1/2}(2 pi xi) / (pi xi)^{p+1/2}
    xi = np.atleast_1d(np.asarray(xi, dtype=float))
    out = np.empty_like(xi); small = np.abs(xi) < 1e-8
    z = np.pi*np.abs(xi[~small])
    out[~small] = np.sqrt(np.pi)*gamma(p+1)*jv(p+0.5, 2*z)/z**(p+0.5)
    out[small] = 2*(wq*Fq).sum()
    return out
F0, Fh0 = 1.0, Fhat(0.0)[0]
Xmax = 400
ns = np.arange(1, Xmax+1)
lhs_atoms = Fhat(ns.astype(float)).sum()
for eps in [0.05+0.0j, 0.03+0.04j, -0.02+0.05j]:
    # D on (m, m+1): 2 Re( eps x^{rho-1} S_m ),  S_m = sum_{n<=m} n^{-rho}
    S = np.cumsum(ns.astype(complex)**(-rho))
    xg2, wg2 = leggauss(40)
    tot = 0.0
    for m in range(1, Xmax):
        x = m + 0.5*(xg2+1); w = 0.5*wg2
        D = 2*np.real(eps * x**(rho-1) * S[m-1])
        tot += (Fhat(x)*D*w).sum()
    lhs = lhs_atoms + tot
    r = Phi(1.0, eps).real; B1 = (Phi(0.0, eps)/Phi(1.0, eps)).real
    rhs = 0.5*r*(F0 - B1*Fh0)
    print("eps=%s  LHS=%.10f  RHS=%.10f  diff=%.2e   (B(1)=%.6f, r=%.6f)" % (eps, lhs, rhs, lhs-rhs, B1, r))
