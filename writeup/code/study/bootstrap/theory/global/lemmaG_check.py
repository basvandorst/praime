# Lemma G identity on the non-self-dual example zeta_P = zeta * (1 + eps/(s-rho) + c.c.):
#   <tau_P, G2> = 0  with  G2(xi) = sin^2(pi xi) / (pi^2 (xi^2 - 1)),
# i.e. int_1^inf G2(x) D(x) dx = 0 (the atoms at integers are killed by G2).  Positivity of tau_P fails,
# as Lemma G predicts, because D changes sign.  Also checks FT(F2) = G2 for F2 = -sin(2 pi |x|)/(2 pi) on [-1,1].
import numpy as np, mpmath as mp
from numpy.polynomial.legendre import leggauss
rho = complex(0.5, float(mp.im(mp.zetazero(1))))
X = 3000
ns = np.arange(1, X+1)
xg, wg = leggauss(60)
def G2(x): return np.sin(np.pi*x)**2/(np.pi**2*(x**2-1))
for eps in [0.05+0.0j, 0.03+0.04j]:
    S = np.cumsum(ns.astype(complex)**(-rho))
    tot = 0.0; Dmin = 1e9
    for m in range(1, X):
        x = m + 0.5*(xg+1); w = 0.5*wg
        D = 2*np.real(eps * x**(rho-1) * S[m-1])
        Dmin = min(Dmin, D.min())
        tot += (G2(x)*D*w).sum()
    cinf = 2*np.real(eps/(1-rho))
    tail = cinf/(2*np.pi**2*X)
    print("eps=%s  int_1^X G2 D = %.3e, + tail %.3e = %.3e   (D_inf=%.4f, min D = %.4f)" % (eps, tot, tail, tot+tail, cinf, Dmin))
# FT check of F2
xg2, wg2 = leggauss(2000); x = 0.5*(xg2+1); w = 0.5*wg2
for xi in [0.0, 0.3, 1.0, 1.5, 2.7]:
    ft = 2*np.sum(w*(-np.sin(2*np.pi*x)/(2*np.pi))*np.cos(2*np.pi*x*xi))
    g = G2(xi) if abs(xi-1) > 1e-9 else 0.0
    print("xi=%.1f  FT(F2)=%.12f  G2=%.12f" % (xi, ft, g))
