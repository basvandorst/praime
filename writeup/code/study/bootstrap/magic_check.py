# Independent numerical test of the magic function printed in Johns, Part II, Theorem 4 (Zenodo, Oct 2026):
#   G(x) = C sin^2(pi x) sqrt(x) sum_n a_n K_0(4 pi sqrt(n x))  (x > 1),  G(1) = C/(32 pi^2),
#   a_n = n S_n + delta_{n1}/4,  S_n = sum_c (1/c)[S(1,n;c) dJ_2(4 pi sqrt n / c) - S(-1,n;c) dI_2(4 pi sqrt n / c)],
#   g(u) = sum_m d(m) m^{-1/2} G(m e^u)   (u >= 0, even).
# Tests with OUR functional W(g) = 2 int g cosh(u/2) + A(g): is W(g) = 0, is the cosine transform >= 0 and zero at zeta zeros?
import numpy as np, math
from scipy import special as sp
from mpmath import mp, mpf
def dord(f, z, h=0.02):                      # derivative with respect to the order at 2, 4-point stencil
    return (-f(2+2*h, z) + 8*f(2+h, z) - 8*f(2-h, z) + f(2-2*h, z))/(12*h)
def kloos(m, n, c):
    if c == 1: return 1.0
    d = np.array([x for x in range(1, c) if math.gcd(x, c) == 1]); dinv = np.array([pow(int(x), -1, c) for x in d])
    return float(np.sum(np.cos(2*np.pi*((m*d + n*dinv) % c)/c)))
def a_small(n, cmax):
    z = 4*np.pi*np.sqrt(n); s = 0.0
    for c in range(1, cmax+1):
        w = z/c
        s += (kloos(1, n, c)*dord(sp.jv, w) - kloos(-1, n, c)*dord(sp.iv, w))/c
    return n*s + (0.25 if n == 1 else 0.0)
NS = 40
a_s = np.array([a_small(n, 400 if n <= 5 else 60) for n in range(1, NS+1)])
print('a_1..a_5 (ours)  :', ' '.join('%.6f' % v for v in a_s[:5]))
print('a_1..a_5 (paper) : 4581.370553 1025727.23492 62590223.68999 1968035704.4549 40671323483.6182')
mp.dps = 30; z1 = 4*mp.pi
print('mpmath c=1 term of a_1:', mp.nstr(mp.diff(lambda m: mp.besselj(m, z1), 2) - mp.diff(lambda m: mp.besseli(m, z1), 2), 15), ' stencil:', dord(sp.jv, float(z1)) - dord(sp.iv, float(z1)))
C = 880.6617
NB = 1500000
nb = np.arange(1, NB+1, dtype=float); zb = 4*np.pi*np.sqrt(nb)
scaled = -nb*dord(sp.ive, zb)                 # a_n e^{-z_n}, c = 1 term (J part and c >= 2 negligible for n > NS)
scaled[:NS] = a_s*np.exp(-zb[:NS])
def Gser(y):                                  # sqrt(y) sum_n a_n K_0(4 pi sqrt(n y)), y > 1
    N = int(min(NB, max(60, (40/(4*np.pi*(np.sqrt(y)-1)))**2 + 10)))
    z = zb[:N]
    return np.sqrt(y)*float(np.sum(scaled[:N]*sp.kve(0, z*np.sqrt(y))*np.exp(-z*(np.sqrt(y)-1))))
dm = [0, 1, 2, 2, 3, 2, 4, 2, 4, 3, 4, 2, 6, 2, 4, 4, 5]
def g(u):
    if u == 0: return C/(32*np.pi**2)
    x = np.exp(u); s = 0.0
    for m in range(1, 17):
        y = m*x
        if 4*np.pi*np.sqrt(y) > 60 and m > 1: break
        s += dm[m]/np.sqrt(m)*C*np.sin(np.pi*y)**2*Gser(y)
    return s
g0 = g(0)
print('g(0) = C/(32 pi^2) =', g0)
print('g(log x)/g(0), paper function vs my LP dual (support log 7.5):')
for x, lp in [(1.02, .9964), (1.1, .9182), (1.201, .7280), (1.327, .4603), (1.502, .1821), (1.748, .02266), (2.504, 1.586e-3), (3.502, 4.499e-5), (4.503, 2.488e-6), (5.498, 4.72e-7), (6.0, 2.26e-7)]:
    print('   x=%.3f  %.4e   %.4e' % (x, g(np.log(x))/g0, lp))
# quadrature: small-u part by a polynomial fit through g(0) and samples on [0.01, 0.08]
us = np.linspace(0.01, 0.08, 15); gs = np.array([g(u) for u in us])
for deg in (5, 6, 7):
    V = np.vander(us, deg+1, increasing=True)[:, 1:]
    co = np.linalg.lstsq(V, gs - g0, rcond=None)[0]; print('   fit degree', deg, 'max residual', np.abs(V@co - (gs-g0)).max(), 'slope at 0:', co[0]/g0)
V = np.vander(us, 7, increasing=True)[:, 1:]; co = np.linalg.lstsq(V, gs - g0, rcond=None)[0]
gfit = lambda u: g0 + sum(c*u**(k+1) for k, c in enumerate(co))
U0 = 0.04
def panel(f, a, b, n=40):
    xg, wg = np.polynomial.legendre.leggauss(n); u = (b-a)/2*xg + (a+b)/2
    return (b-a)/2*float(np.sum(wg*np.array([f(v) for v in u])))
edges = [U0, 0.1, 0.2, 0.35, 0.5, 0.693147, 0.9, 1.098612, 1.25, 1.386294, 1.5, 1.609438, 1.8, 1.94591, 2.1, 2.3, 2.6]
cache = {}
def gg(u):
    if u not in cache: cache[u] = gfit(u) if u < U0 else g(u)
    return cache[u]
def integ(f):
    return panel(f, 0, U0, 30) + sum(panel(f, a, b, 30) for a, b in zip(edges[:-1], edges[1:]))
pole = 4*integ(lambda u: gg(u)*np.cosh(u/2))
arch = -(np.euler_gamma + np.log(np.pi))*g0 + integ(lambda u: (2*np.exp(-2*u)*g0 - 2*np.exp(-u/2)*gg(u))/(1-np.exp(-2*u))) - g0*np.log(1-np.exp(-2*edges[-1]))
print('pole term / g(0) = %.9f   archimedean term / g(0) = %.9f   W(g)/g(0) = %.3e' % (pole/g0, arch/g0, (pole+arch)/g0))
ghat = lambda t: 2*integ(lambda u: gg(u)*np.cos(t*u))
h0 = ghat(0.0)
print('cosine transform, relative to its value at 0:')
ts = np.linspace(0, 60, 601); vals = np.array([ghat(t) for t in ts])/h0
print('   min over t in [0, 60] (step 0.1): %.3e at t = %.1f' % (vals.min(), ts[vals.argmin()]))
for gam in [14.134725142, 21.022039639, 25.010857580, 30.424876126, 32.935061588, 37.586178159, 40.918719012, 43.327073281, 48.005150881, 49.773832478]:
    print('   zero %.6f : %.3e    midpoints -1, +1: %.3e %.3e' % (gam, ghat(gam)/h0, ghat(gam-1)/h0, ghat(gam+1)/h0))

# Truncation to compact support: g_L = g * w_L with w_L = (v*v)/(v*v)(0), v(u) = (1 - (2u/L)^2)_+ .
# w_L >= 0 and its transform is >= 0, so g_L is again doubly positive, with support [-L, L]:
# W(g_L)/g_L(0) is an upper bound (numerical) for the convex optimum at support L.
def window(L):
    b = L/2; xg, wg = np.polynomial.legendre.leggauss(6)
    def vv(u):
        u = abs(u)
        if u >= 2*b: return 0.0
        lo, hi = u - b, b; s = (hi-lo)/2*xg + (hi+lo)/2
        return (hi-lo)/2*float(np.sum(wg*(1 - s**2/b**2)*(1 - (u-s)**2/b**2)))
    v0 = vv(0.0)
    return lambda u: vv(u)/v0
print('truncated magic function, support L = log X:  X, W(g_L)/g_L(0)')
for X in [3, 4, 5, 8, 20, 100, 1e4, 1e8]:
    w = window(np.log(X))
    pole = 4*integ(lambda u: gg(u)*w(u)*np.cosh(u/2))
    arch = -(np.euler_gamma + np.log(np.pi))*g0 + integ(lambda u: (2*np.exp(-2*u)*g0 - 2*np.exp(-u/2)*gg(u)*w(u))/(1-np.exp(-2*u))) - g0*np.log(1-np.exp(-2*edges[-1]))
    print('   X=%-8g  %.4e' % (X, (pole+arch)/g0))
