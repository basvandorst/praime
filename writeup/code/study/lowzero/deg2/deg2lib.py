"""Two-point function of the zeros of the L-function of an elliptic curve E/Q (ratios recipe; Conrey-Snaith method as
written for the Selberg class in D. J. Smith's thesis, arXiv:1607.03297, Heuristic 1) and the number variance Sigma2(l)
in unfolded units. Self-contained: numpy, scipy, mpmath, cypari2.

Normalisation: L(s,E) = sum lambda(n) n^-s, lambda(n) = a_n/sqrt(n); Lambda(s) = N^{s/2} Gamma_C(s+1/2) L(s,E) =
w Lambda(1-s); zeros 1/2 + i gamma. Local roots: good p: alpha, beta = e^{+-i theta}, alpha+beta = a_p/sqrt(p);
multiplicative p: alpha = a_p/sqrt(p) = +-p^{-1/2}, beta = 0; additive: none (not used: the sample is semistable).
Mean density dbar(t) = L/(2 pi), L = ln(N t^2/(2 pi)^2) (+ O(t^-2)); r = dbar*eps, eps = 2 pi r/L.

  R(r) - 1 = (2/L^2) Re[ J(i eps) + e^{-i eps L} K(i eps) ]                                   (law)
  J(x) = sum_p sum_k (ln p)^2 W_p(k) p^{-k(1+x)},  W_p(k) = (alpha^k + beta^k)^2              (diagonal)
  K(x) = prod_p O_p(x)/(1-1/p)^2,  O_p(x) = sum_{m+h=n+k} a(p^m)a(p^n)mu(p^h)mu(p^k) p^{-(m+h)} p^{x(m+n-m-h)}
       = zeta(1+x) zeta(1-x) L(1+x,Sym2) L(1-x,Sym2)/L(1,Sym2)^2 * Aot(x)                     (off-diagonal)
  J(x) = (zeta'/zeta)'(1+x) + (L'/L)'(1+x,Sym2) + Dcorr(x)
Aot and Dcorr are absolutely convergent Euler products/sums (terms O(p^-2)); L(s,Sym2) from PARI (lfunsympow, analytic
s = PARI s - 1; bad factor at p || N: (1 - p^{-1-s})^{-1}).
R(r) - 1 = -(sin pi r/pi r)^2 + C(r);  Sigma2(l) = Sigma2_GUE(l) + 2 int_0^l (l - r) C(r) dr.
"""
import math
import numpy as np
import mpmath as mp
from scipy.interpolate import CubicSpline

TWO_PI = 2 * math.pi
EMAX = 8.0          # eps grid [0, EMAX]
HGRID = 0.01


def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]


def cexpm1(z):
    z = np.asarray(z, complex); x, y = z.real, z.imag
    return (np.expm1(x) * np.cos(y) - 2 * np.sin(y / 2) ** 2) + 1j * np.exp(x) * np.sin(y)


def clog1p(w):
    w = np.asarray(w, complex)
    return 0.5 * np.log1p(2 * w.real + np.abs(w) ** 2) + 1j * np.arctan2(w.imag, 1 + w.real)


def eps_grid(h=HGRID, emax=EMAX):
    return np.arange(0, emax + h / 2, h)


# ---------------------------------------------------------------- zeta part (curve-independent)
class ZetaGrid:
    """Z2 = eps^2|zeta(1+i eps)|^2, Z2m1 = (Z2-1)/eps^2, regd = Re(zeta'/zeta)'(1+i eps) + 1/eps^2; and for the
    degree-1-like alternative: ln Phi_off(eps) = sum_p log(1 - (p^{i eps}-1)^2/(p-1)^2), d2lnPhiDiag."""

    def __init__(self, h=0.0025, emax=EMAX, pmax=10 ** 6, dps=30):
        mp.mp.dps = dps
        eps = np.arange(0, emax + h / 2, h); ee = eps.copy(); ee[0] = 1e-7
        Z2 = np.empty(len(eps)); Z2m1 = np.empty(len(eps)); regd = np.empty(len(eps))
        for i, e in enumerate(ee):
            s = mp.mpc(1, e); z = mp.zeta(s); z1 = mp.zeta(s, derivative=1); z2 = mp.zeta(s, derivative=2); E = mp.mpf(e)
            a = abs(1j * E * z) ** 2; Z2[i] = float(a); Z2m1[i] = float((a - 1) / E ** 2)
            regd[i] = float(mp.re(z2 / z - (z1 / z) ** 2) + 1 / E ** 2)
        P = primes_upto(pmax).astype(float); lnP = np.log(P)
        lnPhi = np.empty(len(eps), complex); d2 = np.empty(len(eps))
        mm = []
        for m in range(2, 64):
            sel = P ** m < 1e18
            if not sel.any(): break
            mm.append((m, lnP[sel], (m - 1) * lnP[sel] ** 2 * P[sel] ** (-float(m))))
        for j in range(0, len(eps), 40):
            e = ee[j:j + 40][:, None]
            w = -4 * np.sin(e * lnP / 2) ** 2 * np.exp(1j * e * lnP) / (P - 1) ** 2
            lnPhi[j:j + 40] = clog1p(-w).sum(axis=1)
            acc = np.zeros(e.shape[0])
            for m, lp, c in mm:
                acc += (c * np.cos(m * e * lp)).sum(axis=1)
            d2[j:j + 40] = 2 * acc
        Phi = np.exp(lnPhi); Phim1 = cexpm1(lnPhi) / ee ** 2
        self.eps = eps
        self.raw = dict(Z2=Z2, Z2m1=Z2m1, regd=regd, Phi=Phi, Phim1=Phim1, d2lnPhiDiag=d2)
        self.sp = {k: ((CubicSpline(eps, v.real), CubicSpline(eps, v.imag)) if np.iscomplexobj(v) else CubicSpline(eps, v))
                   for k, v in self.raw.items()}

    def __call__(self, key, e):
        s = self.sp[key]
        return s[0](e) + 1j * s[1](e) if isinstance(s, tuple) else s(e)


# ---------------------------------------------------------------- local factors
def local_good(p, lam, x):
    """For good primes p (scalar or column array) with lambda = a_p/sqrt(p) and x = i eps (row array): returns
    Dcorr_p(x) = J_p - (zeta'/zeta)'_p - (L'/L)'_{Sym2,p}  = -(ln p)^2 [q(a2 z) + 2 q(z) + q(b2 z)], q(u) = u^2/(1-u)^2,
    lnratio_p(x) = ln[ O_p(x)/(1-1/p)^2 ] - ln[zeta_p(1+x) zeta_p(1-x) LS_p(1+x) LS_p(1-x)/LS_p(1)^2] (Aot factor),
    and the diagonal pieces J_p,k for k = 1, 2 and >= 3 (for decompositions)."""
    p = np.asarray(p, float); lam = np.asarray(lam, float)
    lp = np.log(p); th = np.arccos(np.clip(lam / 2, -1.0, 1.0))
    a2 = np.exp(2j * th); b2 = np.conj(a2)
    z = np.exp(-(1 + x) * lp)
    q = lambda u: u * u / (1 - u) ** 2
    Dcorr = -(lp ** 2) * (q(a2 * z) + 2 * q(z) + q(b2 * z))
    # O_p(x) by its power series: O_p = sum_j p^{-j(1+x)} c_j(p^x)^2, c_j(y) = A_j y^j - lam A_{j-1} y^{j-1} + A_{j-2} y^{j-2}
    y = np.exp(x * lp)
    A = [np.ones_like(lam), lam]                     # A_j = U_j(lam/2): a(p^j) normalised
    O = np.ones(np.broadcast(p, x).shape, complex); j = 1; pmin = float(np.min(p))
    Y = [np.ones_like(y), y]; zj = z                 # Y[j] = y^j, zj = z^j (incremental powers)
    while True:
        if j + 1 >= len(A): A.append(lam * A[-1] - A[-2])
        if j >= len(Y): Y.append(Y[-1] * y)
        c = A[j] * Y[j] - lam * A[j - 1] * Y[j - 1] + (A[j - 2] * Y[j - 2] if j >= 2 else 0)
        O = O + c * c * zj
        if pmin ** (-j) * (j + 2) ** 2 < 1e-19: break
        j += 1; zj = zj * z
        if len(Y) > 3: Y[j - 3] = None
    lnT = np.log(O) - 2 * np.log1p(-1 / p)
    zp = z; zm = np.exp(-(1 - x) * lp); z0 = 1 / p
    lnden = -(np.log(1 - zp) + np.log(1 - zm)) \
            - (np.log(1 - a2 * zp) + np.log(1 - zp) + np.log(1 - b2 * zp)) \
            - (np.log(1 - a2 * zm) + np.log(1 - zm) + np.log(1 - b2 * zm)) \
            + 2 * (np.log(1 - a2 * z0) + np.log(1 - z0) + np.log(1 - b2 * z0)).real
    W1 = lam * lam; W2 = (lam * lam - 2) ** 2
    J1 = lp ** 2 * W1 * z; J2 = lp ** 2 * W2 * z * z
    Jall = lp ** 2 * (a2 * z / (1 - a2 * z) + 2 * z / (1 - z) + b2 * z / (1 - b2 * z))
    return Dcorr, lnT - lnden, J1, J2, Jall - J1 - J2


def local_mult(p, x):
    """multiplicative p (alpha = +-p^{-1/2}, beta = 0): the same quantities."""
    lp = math.log(p); z = np.exp(-(1 + x) * lp)
    q = lambda u: u * u / (1 - u) ** 2; h = lambda u: u / (1 - u) ** 2
    Dcorr = (lp ** 2) * (-q(z / p) - h(z))
    y = np.exp(x * lp)
    O = 1 + (1 - 1 / y) ** 2 * y ** 2 * np.exp(-(1 + x) * lp) / p / (1 - y ** 2 * np.exp(-(1 + x) * lp) / p)
    lnT = np.log(O) - 2 * math.log1p(-1 / p)
    zp = np.exp(-(1 + x) * lp); zm = np.exp(-(1 - x) * lp)
    lnden = -(np.log(1 - zp) + np.log(1 - zm)) - (np.log(1 - zp / p) + np.log(1 - zm / p)) + 2 * math.log1p(-1 / p ** 2)
    J1 = lp ** 2 * z / p; Jall = lp ** 2 * (z / p) / (1 - z / p)
    return Dcorr, lnT - lnden, J1, 0 * J1, Jall - J1


# ---------------------------------------------------------------- one curve
class Curve:
    """All curve-specific regular functions on the eps grid. ainv: minimal model; needs cypari2."""

    def __init__(self, pari, ainv, pmax=10 ** 5, h=HGRID, emax=EMAX, split=(2, 3, 5, 7)):
        E = pari.ellinit(ainv); self.ainv = list(ainv)
        self.N = int(pari.ellglobalred(E)[0]); fa = pari.factor(self.N)
        self.bad = [int(p) for p in fa[0]]; assert all(int(e) == 1 for e in fa[1]), 'semistable only'
        an = pari.ellan(E, pmax); P = primes_upto(pmax)
        self.ap = {int(p): int(an[int(p) - 1]) for p in P}
        eps = eps_grid(h, emax); self.eps = eps; x = 1j * eps
        x[0] = 1e-6j                                     # eps -> 0 limit taken below by extrapolation
        # Sym2 from PARI (analytic s = 1 + i eps  <->  PARI s = 2 + i eps)
        S = pari.lfunsympow(E, 2); Li = pari.lfuninit(S, [2, 0.05, emax + 0.5], 2)
        lnLS = np.empty(len(eps), complex); dd = np.empty(len(eps), complex)
        for i, e in enumerate(eps):
            s = pari('2+%.12f*I' % e)
            v0 = complex(pari.lfun(Li, s)); v1 = complex(pari.lfun(Li, s, 1)); v2 = complex(pari.lfun(Li, s, 2))
            lnLS[i] = np.log(v0); dd[i] = v2 / v0 - (v1 / v0) ** 2
        self.LS1 = float(np.exp(lnLS[0].real))
        # Euler-product parts
        Dc = np.zeros(len(eps), complex); lnA = np.zeros(len(eps), complex)
        self.Jsplit = {}                                 # per split prime: (J1, J2, J3+) on the grid (complex)
        good = np.array([int(p) for p in P if int(p) not in self.bad])
        lam = np.array([self.ap[int(p)] / math.sqrt(p) for p in good])
        xr = x[None, :]
        cuts = [0] + [int(np.searchsorted(good, b)) for b in (50, 1000)] + list(range(int(np.searchsorted(good, 1000)) + 400, len(good), 400)) + [len(good)]
        for k0, k1 in zip(cuts[:-1], cuts[1:]):              # vectorised over chunks of good primes (small primes apart)
            if k1 <= k0: continue
            d, la, _, _, _ = local_good(good[k0:k1, None], lam[k0:k1, None], xr)
            Dc += d.sum(0); lnA += la.sum(0)
        for p in self.bad:
            d, la, _, _, _ = local_mult(p, x); Dc += d; lnA += la
        for p in split:
            if p in self.bad: self.Jsplit[p] = local_mult(p, x)[2:]
            else: self.Jsplit[p] = tuple(v for v in local_good(p, self.ap[p] / math.sqrt(p), x)[2:])
        lnSA = 2 * (lnLS.real - lnLS[0].real) + lnA        # ln[ |LS(1+i eps)|^2/LS(1)^2 * Aot(i eps) ]
        lnSA[0] = 0.0
        KSA = cexpm1(lnSA) / np.where(eps > 0, eps, 1.0) ** 2   # (|LS|^2/LS(1)^2 Aot - 1)/eps^2, regular at 0
        self.raw = dict(dd=dd, Dc=Dc, lnSA=lnSA, KSA=KSA)
        # eps = 0 values of dd, Dc, KSA by cubic extrapolation from the next 4 grid points (smooth functions)
        for k in ('dd', 'Dc', 'KSA'):
            v = self.raw[k]
            c = np.polyfit(eps[1:5], v[1:5].real, 3); ci = np.polyfit(eps[1:5], v[1:5].imag, 3)
            v[0] = np.polyval(c, 0) + 1j * np.polyval(ci, 0)
        self.sp = {k: (CubicSpline(eps, v.real), CubicSpline(eps, v.imag)) for k, v in self.raw.items()}

    def __call__(self, key, e):
        s = self.sp[key]; return s[0](e) + 1j * s[1](e)


# ---------------------------------------------------------------- quadrature and Sigma2
W = 0.025; LMAX = 10.0; NG = 8
_x, _w = np.polynomial.legendre.leggauss(NG)
NPC = int(round(LMAX / W))
R_NODES = (np.arange(NPC)[:, None] * W + (W / 2) * (_x[None, :] + 1)).ravel()
R_W = np.tile((W / 2) * _w, NPC)
LGRID = np.arange(1, NPC + 1) * W


def cum_sigma(f):
    f = f.reshape(NPC, NG); wr = R_W.reshape(NPC, NG); r = R_NODES.reshape(NPC, NG)
    I0 = np.cumsum((f * wr).sum(1)); I1 = np.cumsum((f * wr * r).sum(1))
    return 2 * (LGRID * I0 - I1)


def sigma2_gue():
    return LGRID + cum_sigma(-np.sinc(R_NODES) ** 2)


def C_full(Z, c, L, r=R_NODES, parts=False, kmode='full'):
    """C(r) of the law for curve c at fixed L. kmode: 'full'; 'k1' (traversals k >= 2 of every good prime carry the
    Sato-Tate mean weight 2 instead of W_p(k)); 'offzeta' (off-diagonal amplitude |zeta|^2 Phi_off prod_{p|N} F_p, the
    degree-1 form, diagonal unchanged)."""
    e = TWO_PI * r / L
    Z2 = Z('Z2', e); Z2m1 = Z('Z2m1', e)
    diag = Z('regd', e) + c('dd', e).real + c('Dc', e).real
    lnSA = c('lnSA', e)
    Kreg = Z2m1 * np.exp(lnSA) + c('KSA', e)               # K(i eps) - 1/eps^2
    if kmode == 'offzeta':
        Phi = Z('Phi', e); Phim1 = Z('Phim1', e); Psi = np.ones_like(e, complex); Psim1 = np.zeros_like(e, complex)
        for p in c.bad:
            fp = f_p(p, e); F = 1 + e ** 2 * fp; Psim1 = Psim1 * F + fp; Psi = Psi * F
        Kreg = np.conj(Z2m1 * Phi * Psi + Phim1 * Psi + Psim1)       # BK's convention pairs it with e^{+2 pi i r}
    off = (np.exp(-1j * TWO_PI * r) * Kreg).real
    out = (2 / L ** 2) * (diag + off)
    if kmode == 'k1':
        out = out - (2 / L ** 2) * k1_shift(c, e)
    if parts:
        return out, (2 / L ** 2) * diag, (2 / L ** 2) * off
    return out


def k1_shift(c, e):
    """Re of sum over good p of (ln p)^2 sum_{k>=2} (W_p(k) - 2) p^{-k(1+i eps)} (exact, primes < 10^4)."""
    out = np.zeros_like(e)
    for p, a in c.ap.items():
        if p > 10 ** 4: break
        if p in c.bad: continue
        lam = a / math.sqrt(p); th = math.acos(max(-1.0, min(1.0, lam / 2))); lp = math.log(p)
        k = 2
        while p ** (-k) * k > 1e-18:
            Wk = (2 * math.cos(k * th)) ** 2
            out += lp ** 2 * (Wk - 2) * p ** (-k) * np.cos(k * e * lp); k += 1
    return out


def f_p(p, e):
    lp = math.log(p); sinc = np.sinc(e * lp / (2 * math.pi))
    return (lp ** 2) * sinc ** 2 / (p - 2 + np.exp(1j * e * lp))


def C_deg1(Z, bad, L, r=R_NODES):
    """degree-1-like alternative: every good prime weight 1 at every traversal; bad primes missing orbits (BK 2013 with
    conductor's primes = bad), density L/2 pi."""
    e = TWO_PI * r / L
    Z2m1 = Z('Z2m1', e); Phi = Z('Phi', e); Phim1 = Z('Phim1', e)
    Psi = np.ones_like(e, complex); Psim1 = np.zeros_like(e, complex)
    for p in bad:
        fp = f_p(p, e); F = 1 + e ** 2 * fp; Psim1 = Psim1 * F + fp; Psi = Psi * F
    rego = Z2m1 * Phi * Psi + Phim1 * Psi + Psim1
    d2psi = np.zeros_like(e)
    for p in bad:
        lp = math.log(p); m = 1
        while p ** m < 1e18:
            d2psi += lp ** 2 * np.cos(m * e * lp) / p ** m; m += 1
    d2psi *= 2
    return (2 / L ** 2) * (Z('regd', e) + (np.exp(1j * TWO_PI * r) * rego).real) - (1 / L ** 2) * (Z('d2lnPhiDiag', e) + d2psi)


def L_of_t(N, t):
    """2 pi dbar(t) = ln N - 2 ln 2 pi + 2 Re psi(1 + i t)"""
    return math.log(N) - 2 * math.log(TWO_PI) + 2 * float(mp.re(mp.digamma(mp.mpc(1, t))))


def theta_E(N, t):
    """theta_E(t) = (t/2) ln N - t ln 2 pi + Im ln Gamma(1 + i t); smooth zero count N(t) ~ theta_E(t)/pi - r/2."""
    from scipy.special import loggamma
    t = np.asarray(t, float)
    return t / 2 * math.log(N) - t * math.log(TWO_PI) + np.imag(loggamma(1 + 1j * t))


def L_nodes(N, t1, t2, n=10):
    x, w = np.polynomial.legendre.leggauss(n)
    t = (t1 + t2) / 2 + (t2 - t1) / 2 * x
    Lt = np.array([L_of_t(N, tt) for tt in t]); ww = w * Lt
    return Lt, ww / ww.sum()


def sigma2_curve(Z, c, t1, t2, kmode='full', n=10):
    """height-averaged (zero-density weight) Sigma2(l) on LGRID for curve c, band [t1, t2]"""
    Ls, ws = L_nodes(c.N, t1, t2, n)
    g = sigma2_gue(); return sum(w * (g + cum_sigma(C_full(Z, c, L, kmode=kmode))) for L, w in zip(Ls, ws))


def sigma2_deg1(Z, N, bad, t1, t2, n=10):
    Ls, ws = L_nodes(N, t1, t2, n)
    g = sigma2_gue(); return sum(w * (g + cum_sigma(C_deg1(Z, bad, L))) for L, w in zip(Ls, ws))


def dsig_diag_k(p, Wk, k, L, l):
    """Sigma2 contribution of traversal k of prime p with weight Wk: (2/pi^2) Wk sin^2(pi k l ln p/L)/(k^2 p^k)"""
    return (2 / math.pi ** 2) * Wk * np.sin(math.pi * k * l * math.log(p) / L) ** 2 / (k * k * p ** k)
