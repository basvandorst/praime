# Truncation of Part II's g (g(u) = Fhat(u/2pi)/2pi) with the cubic B-spline window k(u/L), support [-2L,2L].
# W(g_L) = W(g) - int_0^inf K(u) g(u) (1-k(u/L)) du, K(u) = 4cosh(u/2) - 2e^{-u/2}/(1-e^{-2u}); W(g) = 0 (numerically -4.6e-10*g(0)).
import numpy as np
d = np.load('/tmp/claude-1000/conv_prover/fhat.npz'); xis = d['xi']; Fh = d['Fh']
u = 2*np.pi*xis; g = Fh/(2*np.pi); g0 = g[0]; du = u[1]-u[0]
uu = u[1:]; gg = g[1:]
K = 4*np.cosh(uu/2) - 2*np.exp(-uu/2)/(1-np.exp(-2*uu))
def kspl(x):
    x = np.abs(x); return np.where(x <= 1, 1-1.5*x**2+0.75*x**3, np.where(x <= 2, (2-x)**3/4, 0.0))
Ig = np.sum(uu**2*np.exp(-uu/2)*gg/(1-np.exp(-2*uu)))*du
print("I_g/g(0) = %.6f ; universal c1 = 7 zeta(3) + pi^3/4 = %.6f" % (Ig/g0, 7*1.2020569031595942 + np.pi**3/4))
print("second moment int_0^inf K u^2 g du / g(0) = %.6e  (RH forces <= 0)" % (np.sum(K*uu**2*gg)*du/g0))
for L in [0.25, 0.5, 1, 2, 4, 8]:
    dW = np.sum(K*gg*(1-kspl(uu/L)))*du
    print("L=%5.2f support [-%.2f,%.2f]: W(g_L)/g(0) = %.6e   bound 3 I_g/(L^2 g0) = %.3e" % (L, 2*L, 2*L, -dW/g0, 3*Ig/L**2/g0))
print("g(u)/g(0) at u = log 2, 1, 2, 3:", [float(np.interp(x, u, g)/g0) for x in [np.log(2), 1.0, 2.0, 3.0]])
