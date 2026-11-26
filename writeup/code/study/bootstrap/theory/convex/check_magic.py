# Independent numerical test of Part II's F = Xi^2 H_raw, built ONLY from the x>1 formula of Thm 4(b)/(5.8):
#   Ghat_raw(x) = sin^2(pi x) sqrt(x) sum_n a_n K0(4 pi sqrt(n x))  (x>1),   Fhat(xi) = sum_m d(m) m^-1/2 Ghat(m e^{2 pi xi})  (xi>=0),
# Fhat even.  Then F(t) = 2 int_0^inf Fhat(xi) cos(2 pi xi t) dxi.  Tests: F >= 0? F(gamma_k) = 0 (double)? H_raw(0) = F(0)/Xi(0)^2?
import numpy as np, mpmath as mp
from math import gcd
from scipy.special import ive, kve
mp.mp.dps = 30
def Jd(w): return mp.diff(lambda nu: mp.besselj(nu, w), 2)
def Id(w): return mp.diff(lambda nu: mp.besseli(nu, w), 2)
def kloost(m, n, c):
    if c == 1: return 1.0
    d = np.array([x for x in range(c) if gcd(x, c) == 1]); dinv = np.array([pow(int(x), -1, c) for x in d])
    return float(np.sum(np.cos(2*np.pi*((m*d + n*dinv) % c)/c)))
NEX = 12
aex = []   # a_n e^{-z_n}, exact (c <= 800)
for n in range(1, NEX+1):
    z = 4*mp.pi*mp.sqrt(n); S = mp.mpf(0)
    for c in range(1, 801):
        S += (kloost(1, n, c)*Jd(z/c) - kloost(-1, n, c)*Id(z/c))/c
    a = n*S + (mp.mpf(1)/4 if n == 1 else 0)
    aex.append(float(a*mp.e**(-z)))
NMAX = 4_000_000
nn = np.arange(1, NMAX+1, dtype=float); zz = 4*np.pi*np.sqrt(nn)
asc = nn*((2/zz)*ive(1, zz) - (2/zz**2)*ive(0, zz))          # c = 1, -I2dot part, scaled by e^{-z}
for n in range(1, NEX+1):
    if n >= 8: print("n=%d exact/approx = %.15f" % (n, aex[n-1]/asc[n-1]))
asc[:NEX] = aex
def Ghat(x):
    """Ghat_raw(x) for x>1."""
    s = np.sqrt(x) - 1.0
    N = int(min(NMAX, (40.0/(4*np.pi*s))**2 + 50))
    y = zz[:N]*np.sqrt(x)
    S = np.sum(asc[:N]*kve(0, y)*np.exp(zz[:N] - y))
    return np.sin(np.pi*x)**2*np.sqrt(x)*S
def dfun(m): return sum(1 for k in range(1, m+1) if m % k == 0)
G1 = 1/(32*np.pi**2)
for x in [1.3, 1.1, 1.03, 1.01, 1.005, 1.002]:
    print("x=%g  Ghat_raw(x) = %.10f   (claimed Ghat(1) = %.10f)" % (x, Ghat(x), G1))
def Fhat(xi):
    x = np.exp(2*np.pi*xi); tot = 0.0
    for m in range(1, 60):
        y = m*x
        if np.sqrt(y) - 1 > 4: break
        tot += dfun(m)*m**-0.5*(G1 if abs(y-1) < 1e-15 else Ghat(y))
    return tot
h = 0.0005; xis = np.arange(0, 0.8+h/2, h)
Fh = np.array([Fhat(x) for x in xis])
np.savez('/tmp/claude-1000/conv_prover/fhat.npz', xi=xis, Fh=Fh)
print("Fhat(0)=%.10f Fhat min on grid = %.3e  Fhat(xi_n): n=2: %.3e n=3: %.3e" % (Fh[0], Fh.min(), Fhat(np.log(2)/(2*np.pi)), Fhat(np.log(3)/(2*np.pi))))
w = np.full(xis.size, h); w[0] = h/2
def F(t): return 2*np.sum(w*Fh*np.cos(2*np.pi*xis*t))
def Xi(t):
    s = mp.mpf('0.5') + 1j*mp.mpf(t); return float(mp.re(0.5*s*(s-1)*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)))
print("H_raw(0) = F(0)/Xi(0)^2 = %.10f   (Part II Prop 6.4: 0.0011355097904)" % (F(0)/Xi(0)**2))
for t in [0, 5, 10, 13, 14.134725141734693, 15, 20, 21.022039638771555, 25, 25.010857580145688, 30, 30.424876125859513, 40]:
    Ft = F(t); X = Xi(t)
    print("t=%9.5f  F=% .6e  Xi^2=% .6e  H=F/Xi^2=% .6e" % (t, Ft, X*X, Ft/(X*X) if X*X > 1e-30 else float('nan')))
ts = np.linspace(0, 45, 901); Fs = np.array([F(t) for t in ts])
print("min F on [0,45] =", Fs.min(), " at t =", ts[Fs.argmin()])
