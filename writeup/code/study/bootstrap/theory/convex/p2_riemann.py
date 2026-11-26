# P2 numerics: Riemann's h -> f^ = Xi, g0^ = Xi^2. W(g0)/g0(0) two ways; theta-side check; truncation bound.
import numpy as np, mpmath as mp
from scipy.special import digamma
d = np.load('/tmp/claude-1000/conv_prover/xi_grid.npz'); t = d['t']; X = d['xi']; h = t[1]-t[0]
X2 = X**2
w = np.full(t.size, h); w[0] = h/2          # trapezoid on [0,inf) for even integrand: int_R = 2*sum
def Omega(t):
    return np.array([float(mp.re(mp.digamma(mp.mpf(1)/4 + 1j*mp.mpf(x)/2))) for x in t]) - np.log(np.pi)
Om = Omega(t)
g00 = 2*np.sum(w*X2)/(2*np.pi)
A = 2*np.sum(w*X2*Om)/(2*np.pi)
pole = 2*0.25   # g^(i/2)+g^(-i/2) = 2 xi(0)^2 = 1/2
W = pole + A
print("g0(0) =", g00, " A =", A, " W =", W, " W/g0(0) =", W/g00)
def g0(u):
    u = np.atleast_1d(u)
    return np.array([2*np.sum(w*X2*np.cos(t*x))/(2*np.pi) for x in u])
# prime side
import sympy
N = 60
lam = {}
for n in range(2, N+1):
    f = sympy.factorint(n)
    if len(f) == 1: lam[n] = np.log(list(f)[0])
ps = sorted(lam)
vals = g0(np.log(ps))
contrib = [2*lam[n]*n**-0.5*v for n, v in zip(ps, vals)]
print("prime side 2 sum Lambda n^-1/2 g0(log n) =", sum(contrib), " ratio =", sum(contrib)/g00)
for n, c, v in zip(ps[:12], contrib[:12], vals[:12]): print("  n=%d g0(log n)/g0(0)=%.6e contrib/g0(0)=%.6e" % (n, v/g00, c/g00))
# theta side: Theta0(y) = sum 2(2 pi^2 n^4 y^4 - 3 pi n^2 y^2) e^{-pi n^2 y^2}; f0(u) = e^{u/2} Theta0(e^u); f0^ should be Xi
def Theta0(y):
    y = np.asarray(y, float); s = np.zeros_like(y)
    for n in range(1, 40):
        a = np.pi*n*n*y*y
        s += 2*(2*a*a - 3*a)*np.exp(-a)
    return s
uu = np.linspace(-4, 4, 16001); du = uu[1]-uu[0]
f0 = np.exp(uu/2)*Theta0(np.exp(uu))
for tt in [0.0, 3.0, 14.134725141734693, 20.0]:
    print("f0^(%g) = %.12e   Xi = %.12e" % (tt, np.sum(f0*np.cos(tt*uu))*du, float(mp.re(0.5*(0.5+1j*tt)*(-0.5+1j*tt)*mp.pi**(-(0.5+1j*tt)/2)*mp.gamma((0.5+1j*tt)/2)*mp.zeta(0.5+1j*tt)))))
print("f0 even? max|f0(u)-f0(-u)| =", np.max(np.abs(f0-f0[::-1])), " min f0 =", f0.min())
# g0 by autocorrelation on the u grid
for x in [0.0, np.log(2), np.log(3)]:
    k = int(round(x/du)); 
    gx = np.sum(f0[k:]*f0[:f0.size-k])*du if k>0 else np.sum(f0*f0)*du
    print("autocorr g0(%.4f) = %.12e  via Xi^2: %.12e" % (x, gx, g0(x)[0]))
# truncation: g_T = g0 * k(u/T), k = cubic B-spline/M(0); W(g_T) <= W + (3/T^2) int_0^inf e^{-u/2} g0 u^2/(1-e^{-2u}) du
uq = np.linspace(1e-6, 12, 6001); dq = uq[1]-uq[0]
gq = g0(uq)
I2 = np.sum(np.exp(-uq/2)*gq*uq**2/(1-np.exp(-2*uq)))*dq
print("min g0 on [0,12] =", gq.min(), "  I2 =", I2, " I2/g0(0)=", I2/g00)
for T in [1, 2, 3, 5, 10]:
    kk = np.where(uq/T <= 1, 1 - 1.5*(uq/T)**2 + 0.75*(uq/T)**3, np.where(uq/T <= 2, (2-uq/T)**3/4, 0.0))
    Kker = 4*np.cosh(uq/2) - 2*np.exp(-uq/2)/(1-np.exp(-2*uq))
    corr = np.sum(Kker*gq*(1-kk))*dq
    print("T=%g support [-%g,%g]: exact W(g_T)/g(0) = %.6f ; bound W/g0(0)+3 I2/(T^2 g0(0)) = %.6f" % (T, 2*T, 2*T, (W-corr)/g00, W/g00 + 3*I2/T**2/g00))
