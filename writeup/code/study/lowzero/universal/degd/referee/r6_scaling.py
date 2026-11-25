"""Referee check 6 (illustration of G2's threshold, d = 3, Sym^2 11a): smooth bump window W on [a, b], H = W * g_1,
at heights scaled by 2: [20,60], [40,120], [80,240]. Compares
  Z = sum over zeros of H(gamma) e^{2 i theta(gamma)}  (PARI zeros, GRH numerically),
  Plead = -A_3 Q^{-1/6} sum_n Lambda(n) n^{-1/2+1/6} e(-3 (n/Q)^{1/3}) H(t_n*)       (n = p, p^2),
  P1 = same with the first-order bracket H(1 + i kappa/t) + (i/6) H' + (i t/6) H''.
Prints |Z - Plead| (actual leading-kernel error), |P1 - Plead| (actual first-order correction) and the triangle bound
B1 = |A_3| Q^{-1/6} sum |Lambda(n)| n^{-1/2+1/6} |bracket - H| (the bound that G2 uses; expected to grow like T^{(d-1)/2} = T)."""
import math, cmath, time
import numpy as np
from scipy.special import loggamma
import cypari2
P = cypari2.Pari(); P.allocatemem(900 * 10**6)
PI = math.pi
Q, mu, d = 121.0, [1.0, 1.0, 2.0], 3
phi0 = PI / 2 * sum(mu) - d * PI / 4
kap = (d * d - 1) / (24 * d) + 0.5 * sum(m - m * m for m in mu)
A = cmath.exp(1j * (phi0 + PI / 4)) / math.sqrt(d)
def lnGR(z): return -(z / 2) * math.log(PI) + loggamma(z / 2)
def e2it(t):
    s = 0.5 + 1j * np.asarray(t, float); v = 1j * np.asarray(t, float) * math.log(Q) + 0j
    for m in mu: v = v + lnGR(s + m) - lnGR(1 - s + m)
    return np.exp(v)
def Hgrid(a, b):
    """H = W * g_1 and derivatives on a grid (step 0.005), W = exp(1 - 1/(1-x^2)) bump on [a,b]."""
    h = 0.005
    u = np.arange(a - 10, b + 10, h)
    x = (2 * u - a - b) / (b - a)
    W = np.where(np.abs(x) < 1, np.exp(1 - 1 / np.maximum(1e-300, 1 - x * x)), 0.0)
    z = np.arange(-8, 8 + h / 2, h); g = np.exp(-z * z) / math.sqrt(PI) * h
    H = np.convolve(W, g, mode='same')
    H1 = np.gradient(H, h); H2 = np.gradient(H1, h)
    return u, H, H1, H2
E = P.ellinit([0, -1, 1, -10, -20]); LS = P.lfunsympow(E, 2)
t0 = time.time()
zs = np.array([float(x) for x in P.lfunzeros(LS, 260)])
print(f'zeros up to 260: {len(zs)} ({time.time()-t0:.1f}s)', flush=True)
xmax = int(Q * (252 / (2 * PI)) ** 3) + 10
t0 = time.time()
pr = np.array([int(x) for x in P(f'primes([2,{xmax}])')], dtype=np.int64)
ap = np.array([int(x) for x in P(f'my(E=ellinit([0,-1,1,-10,-20])); apply(p->ellap(E,p), primes([2,{xmax}]))')], dtype=float)
print(f'{len(pr)} primes up to {xmax}, a_p in {time.time()-t0:.1f}s', flush=True)
x = ap / np.sqrt(pr)
lam = x * x - 1
good = pr != 11
LamP = np.where(good, lam, 1 / 11.0) * np.log(pr)
# all prime powers p^k, k >= 2: sum alpha^k = beta^{2k} + 1 + beta^{-2k} = 1 + 2 T_{2k}(x/2); at 11: 11^-k
extra_n, extra_L = [], []
for i, p in enumerate(pr):
    p = int(p)
    if p * p > xmax: break
    cc = x[i] / 2; Tm = [1.0, cc]
    k = 2; q = p * p
    while q <= xmax:
        while len(Tm) <= 2 * k: Tm.append(2 * cc * Tm[-1] - Tm[-2])
        extra_n.append(float(q)); extra_L.append((1 + 2 * Tm[2 * k]) * math.log(p) if p != 11 else 11.0 ** -k * math.log(11))
        k += 1; q *= p
ns = np.concatenate([pr.astype(float), np.array(extra_n)]); Lam = np.concatenate([LamP, np.array(extra_L)])
print(f'{len(extra_n)} higher prime powers added', flush=True)
cL = 1.9447054
for (a, b) in ((20.0, 60.0), (40.0, 120.0), (80.0, 240.0)):
    u, H, H1, H2 = Hgrid(a, b)
    Hz = np.interp(zs, u, H)
    Z = np.sum(Hz * e2it(zs))  # negative zeros contribute < 1e-100
    ts = 2 * PI * (ns / Q) ** (1 / 3)
    sel = (ts > a - 9) & (ts < b + 9)
    tt, LL, nn = ts[sel], Lam[sel], ns[sel]
    h0, h1, h2 = np.interp(tt, u, H), np.interp(tt, u, H1), np.interp(tt, u, H2)
    ker = LL * nn ** (-0.5 + 1 / 6) * np.exp(-2j * PI * 3 * (nn / Q) ** (1 / 3))
    corr = h0 * 1j * kap / tt + 1j / 6 * h1 + 1j * tt / 6 * h2
    Plead = -A * Q ** (-1 / 6) * np.sum(ker * h0)
    P1 = Plead - A * Q ** (-1 / 6) * np.sum(ker * corr)
    B1 = abs(A) * Q ** (-1 / 6) * np.sum(np.abs(LL) * nn ** (-0.5 + 1 / 6) * np.abs(corr))
    Babs = abs(A) * Q ** (-1 / 6) * np.sum(np.abs(LL) * nn ** (-0.5 + 1 / 6) * np.abs(h0))
    intH = np.sum(H) * (u[1] - u[0])
    print(f'[{a:.0f},{b:.0f}]: zeros {Z:.4f}; law -c_L intH/2pi = {-cL*intH/(2*PI):.4f}; leading {Plead:.4f}; first order {P1:.4f}\n'
          f'    |Z - Plead| = {abs(Z-Plead):.4f}, |Z - P1| = {abs(Z-P1):.4f}, actual |first-order correction| = {abs(P1-Plead):.4f}; '
          f'triangle bound of the correction B1 = {B1:.2f}; sum |leading terms| = {Babs:.1f}; #n = {len(nn)}', flush=True)
