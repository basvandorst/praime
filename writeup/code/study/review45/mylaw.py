# Independent implementation of the Conrey-Snaith 2008 triple correlation (eq. Jzabb, Qprime, Pprime) for zeta,
# adapted to primitive chi mod k by deleting p | k from every prime sum/product, zeta -> L(., chi0), (k/phi)^2 in P1.
import numpy as np
from scipy.interpolate import CubicSpline
from sympy import primerange, factorint

SCR = __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/'
_t = np.load(SCR + 'ztab.npz')
_ys = np.concatenate([-_t['ys'][::-1], _t['ys']])
def _mk(g):
    gg = np.concatenate([np.conj(g[::-1]), g])
    return CubicSpline(_ys, gg.real), CubicSpline(_ys, gg.imag)
_GZ = _mk(_t['gz']); _GL = _mk(_t['gl']); _GD = _mk(_t['gd'])
def _ev(sp, y): return sp[0](y) + 1j*sp[1](y)

PR_ALL = np.array(list(primerange(2, 1000001)), dtype=float)
SMALL = PR_ALL[PR_ALL <= 13]

def zeta1(y):      # zeta(1+iy), y real (nonzero)
    return _ev(_GZ, y) + 1/(1j*y)
def zlogd(y):      # zeta'/zeta(1+iy)
    return _ev(_GL, y) - 1/(1j*y)
def zlogdd(y):     # (zeta'/zeta)'(1+iy)
    return _ev(_GD, y) - 1/y**2

# A over primes 13 < p <= 1e5, tabulated; small primes exact
_PA = PR_ALL[(PR_ALL > 13) & (PR_ALL <= 1e5)]
_yA = np.linspace(-6.6, 6.6, 2641)
_lA = np.zeros(len(_yA), complex)
for i0 in range(0, len(_yA), 200):
    yy = _yA[i0:i0+200, None]
    z = np.exp(-1j*yy*np.log(_PA)[None, :])
    _lA[i0:i0+200] = np.log1p(-(1 - z)**2/(_PA - 1)**2).sum(1)
_SA = (CubicSpline(_yA, _lA.real), CubicSpline(_yA, _lA.imag))
def Afac(p, y):    # local factor of A at p, x = iy
    z = np.exp(-1j*y*np.log(p))
    return 1 - (1 - z)**2/(p - 1)**2
def A_k(y, bad):
    v = np.exp(_ev(_SA, y))
    for p in SMALL:
        if p not in bad: v = v*Afac(p, y)
    for p in bad:
        if p > 13: v = v/Afac(p, y)
    return v

def Lchi0(y, bad):   # L(1+iy, chi0)
    v = zeta1(y)
    for p in bad: v = v*(1 - p**(-1.0)*np.exp(-1j*y*np.log(p)))
    return v
def Llogd(y, bad):   # L'/L(1+iy, chi0)
    v = zlogd(y)
    for p in bad: v = v + np.log(p)/(p*np.exp(1j*y*np.log(p)) - 1)
    return v
def Llogdd(y, bad):  # (L'/L)'(1+iy, chi0)
    v = zlogdd(y)
    for p in bad:
        q = p*np.exp(1j*y*np.log(p))
        v = v - np.log(p)**2*q/(q - 1)**2
    return v

def P1(y, L, bad, kphi2):
    return np.exp(-1j*L*y)*kphi2*Lchi0(y, bad)*Lchi0(-y, bad)*A_k(y, bad)

def B1(yx, yy, bad, P=3000, tail=True):
    """CS 2008 B1(x, y), x = i yx, y = i yy, sum over primes p <= P not dividing k, plus PNT tail."""
    pr = PR_ALL[PR_ALL <= P]
    pr = pr[~np.isin(pr, list(bad))]
    out = np.zeros(yx.shape, complex)
    lp = np.log(pr)
    for i0 in range(0, yx.size, 1500):
        ax = yx.ravel()[i0:i0+1500, None]; ay = yy.ravel()[i0:i0+1500, None]
        ex = np.exp(-1j*ax*lp); ey = np.exp(-1j*ay*lp)   # p^{-x}, p^{-y}
        p = pr[None, :]
        num = (1 - ex)*(1 - ex - ey + ey/p)*lp
        den = (1 - ey/(p*ex))*(1 - ey/p)*(1 - 2/p + ex/p)*p**2*ex/ey
        out.ravel()[i0:i0+1500] = (num/den).sum(1)
    if tail:
        # b1 ~ ln p p^{-2} * sum c_j p^{-i s_j}; sum_{p>P} ln p p^{-2-is} ~ P^{-1-is}/(1+is)
        a = yx; b = yy
        terms = [(1, b - a), (-2, b), (-1, 2*b - a), (1, a + b), (1, 2*b)]
        for c, s in terms:
            out = out + c*np.exp(-(1 + 1j*s)*np.log(P))/(1 + 1j*s)
    return out

def Bsum(y, bad, P=1000000):  # sum_{p not | k} (ln p/(p^{1+x}-1))^2, x = iy (1-D)
    pr = PR_ALL[PR_ALL <= P]
    pr = pr[~np.isin(pr, list(bad))]
    lp = np.log(pr)
    out = np.zeros(y.shape, complex)
    for i0 in range(0, y.size, 50):
        yy = y.ravel()[i0:i0+50, None]
        out.ravel()[i0:i0+50] = ((lp/(pr*np.exp(1j*yy*lp) - 1))**2).sum(1)
    return out

def badset(D):
    return [float(p) for p in factorint(abs(D))]
def kphi2(D):
    v = 1.0
    for p in factorint(abs(D)): v *= p/(p - 1)
    return v**2

def T3_P(t, s, L, bad, kp2, P=3000, gue=False):
    """connected 3-pt function minus Q-terms, unfolded positions (0, t, s); returns real array."""
    c = 2*np.pi/L
    u = np.zeros_like(t); v = c*t; w = c*s
    tot = np.zeros(t.shape, complex)
    for (z, x, yv) in ((w, u, v), (v, w, u), (u, w, v)):
        # J*(-iz; ix, iy) without Q
        X = x - z; Yv = yv - z      # a+b1 = iX, a+b2 = iYv, b2-b1 = i(Yv-X)
        if gue:
            p1X = np.exp(-1j*L*X)*(-1/(1j*X)**2); p1Y = np.exp(-1j*L*Yv)*(-1/(1j*Yv)**2)
            P3a = -1/(1j*Yv) + 1/(1j*(Yv - X))
            P3b = -1/(1j*X) + 1/(1j*(X - Yv))
        else:
            p1X = P1(X, L, bad, kp2); p1Y = P1(Yv, L, bad, kp2)
            P3a = B1(X, Yv, bad, P) + Llogd(Yv, bad) - Llogd(Yv - X, bad)
            P3b = B1(Yv, X, bad, P) + Llogd(X, bad) - Llogd(X - Yv, bad)
        tot = tot - p1X*P3a - p1Y*P3b
    return 2*tot.real/L**3

def kQ_closed(l, L, bad, P=1000000, M=40):
    """kappa[3]_Q = -(6/pi^3) sum_{p not | k} sum_{m,n} sin m th sin n th sin (m+n) th / (m n (m+n) p^(m+n))"""
    pr = PR_ALL[PR_ALL <= P]; pr = pr[~np.isin(pr, list(bad))]
    th = np.pi*l*np.log(pr)/L
    tot = np.zeros_like(pr)
    for m in range(1, M):
        for n in range(1, M - m + 1):
            w = pr**(-float(m + n))
            if w[0] < 1e-18: continue
            tot += np.sin(m*th)*np.sin(n*th)*np.sin((m + n)*th)/(m*n*(m + n))*w
    return -(6/np.pi**3)*tot.sum()

def T2(r, L, bad, kp2):
    y = 2*np.pi*r/L
    P2 = Llogdd(y, bad) - Bsum(y, bad)
    return 2*(P1(y, L, bad, kp2) + P2).real/L**2

# quadrature nodes
def gl(n): return np.polynomial.legendre.leggauss(n)
def nodes(lmax=8.0, piece=0.25, ng=8):
    xg, wg = gl(ng)
    S = []; WS = []
    for a in np.arange(0, lmax, piece):
        S.append(a + piece*(xg + 1)/2); WS.append(piece*wg/2)
    S = np.concatenate(S); WS = np.concatenate(WS)
    T = []; ST = []; WT = []; IDX = []
    for i, s in enumerate(S):
        half = s/2; npc = max(1, int(np.ceil(half/piece)))
        h = half/npc
        for j in range(npc):
            T.append(j*h + h*(xg + 1)/2); WT.append(h*wg/2); IDX.append(np.full(ng, i))
    return S, WS, np.concatenate(T), np.concatenate(WT), np.concatenate(IDX)

def kf3_P(ls, L, bad, kp2, P=3000, gue=False, nd=None):
    S, WS, T, WT, IDX = nd if nd is not None else nodes()
    f = T3_P(T, S[IDX], L, bad, kp2, P, gue)
    F = 2*np.bincount(IDX, weights=f*WT, minlength=len(S))   # F(s) = 2 int_0^{s/2} T3 dt
    return np.array([6*np.sum(WS*(l - S)*F*(S < l)) for l in ls])

def sigma2(ls, L, bad, kp2, ng=16):
    xg, wg = gl(ng)
    out = []
    for l in ls:
        R = []; W = []
        for a in np.arange(0, l, 0.25):
            R.append(a + 0.25*(xg + 1)/2); W.append(0.25*wg/2)
        R = np.concatenate(R); W = np.concatenate(W)
        out.append(l + 2*np.sum(W*(l - R)*T2(R, L, bad, kp2)))
    return np.array(out)

def Qsum(X, Y, bad, P=200000):
    pr = PR_ALL[PR_ALL <= P]; pr = pr[~np.isin(pr, list(bad))]; lp = np.log(pr)
    out = np.zeros(X.shape, complex)
    for i0 in range(0, X.size, 200):
        a = X.ravel()[i0:i0+200, None]; b = Y.ravel()[i0:i0+200, None]
        out.ravel()[i0:i0+200] = -(lp**3/((pr*np.exp(1j*a*lp) - 1)*(pr*np.exp(1j*b*lp) - 1))).sum(1)
    return out
def T3_Q(t, s, L, bad):
    c = 2*np.pi/L
    u = np.zeros_like(t); v = c*t; w = c*s
    tot = np.zeros(t.shape, complex)
    for (z, x, yv) in ((w, u, v), (v, w, u), (u, w, v)):
        tot = tot + Qsum(x - z, yv - z, bad)
    return 2*tot.real/L**3
