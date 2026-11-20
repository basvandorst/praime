"""Smooth numbers: Psi(x, y) = #{n <= x : every prime factor of n is <= y} against
(a) Dickman: x rho(u), u = ln x / ln y;  (b) de Bruijn's second-order term: x rho(u) + (1 - gamma) x rho(u - 1) / ln x;
(c) de Bruijn's Lambda(x, y) = x int rho(u - v) d(floor(y^v)/y^v)
    = x [ sum_{m <= x} rho(u - ln m/ln y)/m - int_1^x rho(u - ln t/ln y) floor(t)/t^2 dt ]   (Gauss-Legendre, 6 nodes per unit interval).
rho: rho = 1 on [0, 1], u rho'(u) = -rho(u - 1), trapezoid with step 1e-4.
Predictions written before the run: Lambda within 0.1 % (Saias), de Bruijn's term roughly halves the error of rho,
rho itself off by 10 % or more at u = 3 to 5.   usage: smooth.py x out.txt"""
import sys, math
import numpy as np
x = int(sys.argv[1]); out = open(sys.argv[2], 'w')
def pr(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')
h = 1e-4; U = np.arange(0, 20 + h / 2, h); rho = np.ones_like(U); k1 = int(round(1 / h))
for i in range(k1 + 1, len(U)):
    rho[i] = (U[i - 1] * rho[i - 1] - h * (rho[i - k1] + rho[i - 1 - k1]) / 2) / U[i] if False else rho[i - 1] - h * (rho[i - k1] / U[i] + rho[i - 1 - k1] / U[i - 1]) / 2
def R(u): u = np.asarray(u, float); return np.where(u < 0, 0.0, np.interp(u, U, rho))
pr('rho(2) = %.6f (1 - ln 2 = %.6f), rho(3) = %.6f, rho(4) = %.7f' % (R(2), 1 - math.log(2), R(3), R(4)))
gpf = np.ones(x + 1, np.int64); gpf[0] = 0
sv = np.ones(x + 1, bool); sv[:2] = False
for i in range(2, int(x ** .5) + 1):
    if sv[i]: sv[i * i::i] = False
P = np.nonzero(sv)[0]
for p in P: gpf[p::p] = p
cnt = np.bincount(gpf[1:], minlength=x + 1); cum = np.cumsum(cnt)        # cum[y] = Psi(x, y)
m = np.arange(1, x + 1, dtype=float); lm = np.log(m)
GL, GW = np.polynomial.legendre.leggauss(6)
def Lam(y, Rf=None):
    Rf = Rf or R
    ly = math.log(y); u = math.log(x) / ly
    s1 = float(np.sum(Rf(u - lm / ly) / m))
    f = lambda t: Rf(u - np.log(t) / ly) / t ** 2
    s2 = 0.0
    for lo in range(1, x, 10 ** 6):                      # Gauss-Legendre with 6 nodes on every [a, a + 1]
        a = np.arange(lo, min(x, lo + 10 ** 6), dtype=float)
        s2 += float(np.sum(a * sum(w * f(a + (g + 1) / 2) for g, w in zip(GL, GW)) / 2))
    return x * (s1 - s2)
pr('check: with rho = 1 Lambda must be x: %.6f x (the plain Simpson rule gave 0.995 x at x = 1000)' % (Lam(10, lambda v: np.where(np.asarray(v) < 0, 0.0, 1.0 + 0 * np.asarray(v))) / x))
g = 0.5772156649015329
pr('x = %d' % x)
pr('      y       u        Psi      x rho   ratio | + de Bruijn term   ratio | Lambda(x,y)   ratio')
ys = sorted(set([10, 30, 100, 300, 1000] + [int(round(x ** (1 / u))) for u in (1.5, 2, 2.5, 3, 3.5, 4, 5, 6)]))
for y in ys:
    if y < 2 or y > x: continue
    u = math.log(x) / math.log(y); ps = int(cum[y]); a = x * float(R(u)); b = a + (1 - g) * x * float(R(u - 1)) / math.log(x); c = Lam(y)
    pr('%7d %7.3f %10d %10.1f %7.4f | %16.1f %7.4f | %11.1f %7.4f' % (y, u, ps, a, a / ps, b, b / ps, c, c / ps))
