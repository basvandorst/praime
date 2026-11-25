# Own implementation of the ratios-recipe two-point law for L(s,E) (Smith 2016 Heuristic 1 / BKS Conj. RCd),
# its number variance, alternative templates and explicit-formula (diagonal, cut at X) number variances.
import numpy as np, math, cypari2
from scipy.interpolate import CubicSpline
from scipy.special import digamma

pari = cypari2.Pari(); pari.allocatemem(600 * 10 ** 6); pari.set_real_precision(30)
ZC = np.load(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/zeta_cache.npz')
EPS = ZC['eps']; Z0 = ZC['Z0']; ZD = ZC['ZD']
LS = np.arange(0.5, 10.01, 0.5)
PMAX = 100000
PR = np.array([int(p) for p in pari.primes(pari.primepi(PMAX))])
LNP = np.log(PR)


def Lt(N, t):
    return math.log(N) - 2 * math.log(2 * math.pi) + 2 * np.real(digamma(1 + 1j * np.asarray(t, float)))


def nodes(lo=20.0, hi=400.0, n=12):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (hi - lo) * x + 0.5 * (hi + lo), 0.5 * (hi - lo) * w


def gue_s2(ls=LS, h=0.0005):
    out = []
    for l in ls:
        r = np.arange(0, l + h / 2, h); f = (l - r) * np.sinc(r) ** 2
        out.append(l - 2 * (h / 3) * (f[0] + f[-1] + 4 * f[1:-1:2].sum() + 2 * f[2:-1:2].sum()))
    return np.array(out)


GUE = gue_s2()


class Curve:
    def __init__(self, ainv, N):
        self.N = N; self.E = pari.ellinit(ainv)
        an = pari.ellan(self.E, PMAX)
        a = np.array([int(x) for x in an], dtype=np.int64)
        self.ap = a[PR - 1]
        self.bad = np.array([N % int(p) == 0 for p in PR])
        self.lam = self.ap / np.sqrt(PR)
        # Sym^2 L-function at analytic 1 + i e  (PARI s = 2 + i e)
        Ls = pari.lfunsympow(self.E, 2)
        Li = pari.lfuninit(Ls, [2, 0, 8.3], 2)
        eg = np.arange(0, 8.3, 0.01); v0 = []; v1 = []; v2 = []
        for e in eg:
            s = pari('2+%.10f*I' % e)
            v0.append(complex(pari.lfun(Li, s, 0))); v1.append(complex(pari.lfun(Li, s, 1))); v2.append(complex(pari.lfun(Li, s, 2)))
        sp = lambda v: CubicSpline(eg, np.real(v))(EPS) + 1j * CubicSpline(eg, np.imag(v))(EPS)
        self.S0 = sp(v0); self.S1 = sp(v1); self.S2 = sp(v2)
        self.LS1 = float(pari.lfun(Ls, 2))
        self.Aot = self._aot(pmax=20000)
        self.Az = self._azeta()

    # ---- local pieces
    def W(self, k, idx=None):
        """W_p(k) = (alpha^k + beta^k)^2 for all p < PMAX (good), p^-k (bad)"""
        lam = self.lam; b0 = np.full_like(lam, 2.0); b1 = lam.copy()
        if k == 1: bk = b1
        else:
            for _ in range(k - 1): b0, b1 = b1, lam * b1 - b0
            bk = b1
        W = bk ** 2
        W[self.bad] = PR[self.bad] ** (-float(k))
        return W

    def _aot(self, pmax=PMAX):
        """A_ot(i e) = prod_p O_p(i e) Lx_p(1)^2/(Lx_p(1+ie) Lx_p(1-ie)), Lx_p = zeta_p * Sym2_p (analytic)"""
        logA = np.zeros(len(EPS), complex)
        npr = int(np.searchsorted(PR, pmax, side='right'))
        for i0 in range(0, npr, 400):
            i1 = min(i0 + 400, npr)
            p = PR[i0:i1].astype(float); lam = self.lam[i0:i1]; bad = self.bad[i0:i1]
            J = 60 if p[0] < 50 else (20 if p[0] < 1000 else 7)
            y = np.exp(1j * np.outer(np.log(p), EPS)); yi = 1 / y
            lm = [np.ones_like(lam), lam.copy()]
            for m in range(2, J + 1): lm.append(np.where(bad, lam * lm[-1], lam * lm[-1] - lm[-2]))
            mu1 = -lam; mu2 = np.where(bad, 0.0, 1.0)
            O = np.ones_like(y); yj = np.ones_like(y); yjm1 = None; yjm2 = None; pyi = np.ones_like(y)
            for j in range(1, J + 1):
                yjm2, yjm1, yj = yjm1, yj, yj * y
                pyi = pyi * (yi / p[:, None])
                c = lm[j][:, None] * yj + (mu1 * lm[j - 1])[:, None] * yjm1
                if j >= 2: c = c + (mu2 * lm[j - 2])[:, None] * yjm2
                O = O + pyi * c * c
            pp = p[:, None]; l2 = (lam ** 2 - 2)[:, None]; bd = bad[:, None]
            def logLx(X):
                return np.where(bd, -np.log(1 - X) - np.log(1 - X / pp), -2 * np.log(1 - X) - np.log(1 - l2 * X + X * X))
            X1 = (1.0 / pp) * np.ones_like(y); Xp = yi / pp; Xm = y / pp
            logA += np.sum(np.log(O) + 2 * logLx(X1) - logLx(Xp) - logLx(Xm), axis=0)
        return np.exp(logA)

    def _azeta(self):
        """degree-1 (Dirichlet, modulus N) off-diagonal arithmetic factor: prod_{p not| N} A_p(ie) prod_{p|N} bad factor"""
        logA = np.zeros(len(EPS), complex)
        for i0 in range(0, len(PR), 2000):
            p = PR[i0:i0 + 2000].astype(float); bad = self.bad[i0:i0 + 2000]
            x = 1j * EPS[None, :]; pp = p[:, None]
            good = np.log((1 - pp ** (-1 - x)) * (1 - 2 / pp + pp ** (-1 - x)) / (1 - 1 / pp) ** 2)
            bd = np.log((1 - pp ** (-1 - x)) * (1 - pp ** (-1 + x)) / (1 - 1 / pp) ** 2)
            logA += np.sum(np.where(bad[:, None], bd, good), axis=0)
        return np.exp(logA)

    # ---- J and K on the eps grid, regularised
    def FJ(self):
        """J_RS(ie) + 1/e^2 with J_RS = (log zeta L(Sym2))''(1+ie)"""
        return ZD + (self.S2 / self.S0 - (self.S1 / self.S0) ** 2) + 1 / EPS ** 2

    def FK(self, mode='full'):
        if mode == 'full':
            K = np.abs(Z0) ** 2 * np.abs(self.S0) ** 2 / self.LS1 ** 2 * self.Aot
        elif mode == 'nosym':        # Sym2 ratio replaced by 1 (A_ot kept)
            K = np.abs(Z0) ** 2 * self.Aot
        elif mode == 'deg1':         # Dirichlet-type off-diagonal
            K = np.abs(Z0) ** 2 * self.Az
        return K - 1 / EPS ** 2

    def diag_corr(self, L, ls=LS, kmax=40):
        """number-variance contribution of J - J_RS (absolutely convergent): (2/pi^2) sum [W - k V] sin^2/(k^2 p^k)"""
        out = np.zeros(len(ls))
        for k in range(1, kmax + 1):
            m = PR.astype(float) ** k < 1e30
            W = self.W(k)
            V = W.copy(); V[self.bad] = 1 + PR[self.bad] ** (-float(k))
            c = (W - k * V)
            if k == 1: c[~self.bad] = 0.0
            sel = m & (np.abs(c) / PR.astype(float) ** k > 1e-18)
            if not sel.any(): continue
            ph = np.pi * k * np.outer(ls, LNP[sel]) / L
            out += (2 / np.pi ** 2) * (np.sin(ph) ** 2 * (c[sel] / (k * k * PR[sel].astype(float) ** k))).sum(1)
        return out


def s2_from_F(FJ, FK, L, ls=LS, h=0.0025):
    """Sigma2(l) - GUE(l) = 2 int_0^l (l - r) C(r) dr, C = (2/L^2) Re[FJ + e^{-2 pi i r} FK], e = 2 pi r/L"""
    e0 = np.concatenate([[0.0], EPS])
    def ext(F):                       # extrapolate to e = 0 quadratically from the first 3 points
        f0 = 3 * F[0] - 3 * F[1] + F[2]
        return np.concatenate([[f0], F])
    sJr = CubicSpline(e0, ext(FJ.real)); sKr = CubicSpline(e0, ext(FK.real)); sKi = CubicSpline(e0, ext(FK.imag))
    rmax = ls[-1]; r = np.arange(0, rmax + h / 2, h); e = 2 * np.pi * r / L
    assert e[-1] <= EPS[-1]
    C = (2 / L ** 2) * (sJr(e) + np.cos(2 * np.pi * r) * sKr(e) + np.sin(2 * np.pi * r) * sKi(e))
    out = []
    for l in ls:
        n = int(round(l / h)); rr = r[:n + 1]; f = (l - rr) * C[:n + 1]
        out.append(2 * (h / 3) * (f[0] + f[-1] + 4 * f[1:-1:2].sum() + 2 * f[2:-1:2].sum()))
    return np.array(out)


def diag_X(cv, L, X, ls=LS, kmax=40, tail_to=None):
    """explicit formula, orthogonal lines: (2/pi^2) sum_{p^k <= X} W sin^2(pi k l ln p/L)/(k^2 p^k);
       for X > PMAX the primes PMAX < p <= X are added with the Sato-Tate/Rankin-Selberg mean weight 1 (k = 1)"""
    out = np.zeros(len(ls))
    for k in range(1, kmax + 1):
        pk = PR.astype(float) ** k; sel = pk <= X
        if not sel.any(): break
        W = cv.W(k)[sel]
        ph = np.pi * k * np.outer(ls, LNP[sel]) / L
        out += (2 / np.pi ** 2) * (np.sin(ph) ** 2 * (W / (k * k * pk[sel]))).sum(1)
    if X > PMAX:
        # int_{PMAX}^{X} sin^2(pi l ln u/L)/(u ln u) du, substitution v = ln u
        v = np.linspace(math.log(PMAX), math.log(X), 4001)
        f = np.sin(np.pi * np.outer(ls, v) / L) ** 2 / v
        out += (2 / np.pi ** 2) * np.trapezoid(f, v, axis=1)
    return out


def diag_corr_w(cv, L, wfun, ls=LS, kmax=40):
    """(2/pi^2) sum_{p,k} [w(p,k) - k V_RS(p,k)] sin^2/(k^2 p^k) for a template weight w, relative to (log zeta)''
       (V_RS = 1 for every p): used for the deg1 template whose FJ is (zeta'/zeta)' only"""
    out = np.zeros(len(ls))
    for k in range(1, kmax + 1):
        w = wfun(k); c = w - k * 1.0
        sel = np.abs(c) / PR.astype(float) ** k > 1e-18
        if not sel.any(): continue
        ph = np.pi * k * np.outer(ls, LNP[sel]) / L
        out += (2 / np.pi ** 2) * (np.sin(ph) ** 2 * (c[sel] / (k * k * PR[sel].astype(float) ** k))).sum(1)
    return out


def single(L, p, k, dW, ls=LS):
    return (2 / np.pi ** 2) * dW * np.sin(np.pi * k * ls * math.log(p) / L) ** 2 / (k * k * p ** k)


def gue_completion(L, X, ls=LS):
    """(2/pi^2)[int_{tX}^1 sin^2(pi l tau)/tau dtau + int_1^inf sin^2(pi l tau)/tau^2 dtau], tX = ln X/L"""
    tX = min(math.log(X) / L, 1.0); out = []
    for l in ls:
        t = np.linspace(tX, 1, 20001); a = np.trapezoid(np.sin(np.pi * l * t) ** 2 / t, t) if tX < 1 else 0.0
        # int_1^inf sin^2(pi l tau)/tau^2 = 1/2 int_1^inf (1 - cos 2 pi l tau)/tau^2 = 1/2 [1 - (cos(2 pi l) - 2 pi l (pi/2 - Si(2 pi l)))]
        from scipy.special import sici
        z = 2 * np.pi * l; si, ci = sici(z)
        b = 0.5 * (1 - (np.cos(z) - z * (np.pi / 2 - si)))
        out.append((2 / np.pi ** 2) * (a + b))
    return np.array(out)
