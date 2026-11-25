"""Bogomolny-Keating (2013, arXiv:1307.6010) two-point function of the zeros of a primitive Dirichlet L-function of
conductor k, turned into the number variance Sigma2(l) in unfolded units. Self-contained (numpy, scipy, mpmath).

Paper (eq. numbers of the arXiv source):
  R2(eps) = dbar^2 + R2diag + R2off,  dbar = ln(kE/2pi)/(2pi) =: L/(2pi)
  R2diag  = -(1/4pi^2) d^2/deps^2 ln[ |zeta(1+i eps)|^2 Phi_diag(eps) Psi(eps,k) ]
  Phi_diag = prod_p exp(2 sum_m (1-m)/(m^2 p^m) cos(m eps ln p)),  Psi = exp(-2 sum_{p|k} sum_m cos(m eps ln p)/(m^2 p^m))
  R2off   = (1/4pi^2) e^{i eps L} |zeta(1+i eps)|^2 Phi_off(eps) Psi_off(eps,k) + c.c.
  Phi_off = prod_p (1 - (p^{i eps}-1)^2/(p-1)^2),  Psi_off = prod_{p|k} F_p,  F_p = (p - p^{-i eps})/(p - 2 + p^{i eps})
Unfolded distance r = dbar*eps, eps = 2 pi r/L.  R(r) - 1 = (R2 - dbar^2)/dbar^2 = -(sin pi r/pi r)^2 + C(r).
Sigma2(l) = l + 2 int_0^l (l-r)(R(r)-1) dr = Sigma2_GUE(l) + 2 int_0^l (l-r) C(r) dr.
"""
import math
import numpy as np
import mpmath as mp
from scipy.interpolate import CubicSpline

TWO_PI = 2 * math.pi


def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]


def clog1p(w):
    """accurate complex log(1+w)"""
    w = np.asarray(w, complex)
    re = 0.5 * np.log1p(2 * w.real + np.abs(w) ** 2)
    im = np.arctan2(w.imag, 1 + w.real)
    return re + 1j * im


def cexpm1(z):
    """accurate complex exp(z)-1"""
    z = np.asarray(z, complex)
    x, y = z.real, z.imag
    return (np.expm1(x) * np.cos(y) - 2 * np.sin(y / 2) ** 2) + 1j * np.exp(x) * np.sin(y)


class Grid:
    """Regular parts of the BK two-point function tabulated on eps in [0, eps_max], cubic splines in between."""

    def __init__(self, eps_max=7.6, h=0.0025, pmax=10 ** 6, dps=30, verbose=True):
        mp.mp.dps = dps
        eps = np.arange(0, eps_max + h / 2, h)
        eps_eval = eps.copy(); eps_eval[0] = 1e-7          # the eps -> 0 limit (all functions regular there)
        self.eps = eps
        P = primes_upto(pmax).astype(float); lnP = np.log(P)
        # zeta part: Z2 = eps^2 |zeta(1+i eps)|^2, Z2m1 = (Z2-1)/eps^2, regd = Re (zeta'/zeta)'(1+i eps) + 1/eps^2
        Z2 = np.empty(len(eps)); Z2m1 = np.empty(len(eps)); regd = np.empty(len(eps))
        for i, e in enumerate(eps_eval):
            s = mp.mpc(1, e); z = mp.zeta(s); z1 = mp.zeta(s, derivative=1); z2 = mp.zeta(s, derivative=2)
            ee = mp.mpf(e)
            a = abs(1j * ee * z) ** 2
            Z2[i] = float(a); Z2m1[i] = float((a - 1) / ee ** 2)
            regd[i] = float(mp.re(z2 / z - (z1 / z) ** 2) + 1 / ee ** 2)
        # Phi_off: ln Phi_off = sum_p log(1 - w_p), w_p = (p^{i e}-1)^2/(p-1)^2 = -4 sin^2(e ln p/2) p^{i e}/(p-1)^2
        lnPhi = np.empty(len(eps), complex); d2lnPhiDiag = np.empty(len(eps))
        # Phi_diag: d^2/de^2 ln Phi_diag = 2 sum_p sum_{m>=2} (m-1) (ln p)^2 p^{-m} cos(m e ln p)
        mm = [];
        for m in range(2, 64):
            sel = P ** m < 1e18
            if not sel.any(): break
            mm.append((m, lnP[sel], ((m - 1) * lnP[sel] ** 2 * P[sel] ** (-float(m)))))
        CH = 40
        for j in range(0, len(eps), CH):
            e = eps_eval[j:j + CH][:, None]
            w = -4 * np.sin(e * lnP / 2) ** 2 * np.exp(1j * e * lnP) / (P - 1) ** 2
            lnPhi[j:j + CH] = clog1p(-w).sum(axis=1)
            acc = np.zeros(e.shape[0])
            for m, lp, c in mm:
                acc += (c * np.cos(m * e * lp)).sum(axis=1)
            d2lnPhiDiag[j:j + CH] = 2 * acc
        Phi = np.exp(lnPhi); Phim1 = cexpm1(lnPhi) / eps_eval ** 2
        self.raw = dict(Z2=Z2, Z2m1=Z2m1, regd=regd, Phi=Phi, Phim1=Phim1, d2lnPhiDiag=d2lnPhiDiag)
        self.sp = {}
        for k, v in self.raw.items():
            if np.iscomplexobj(v):
                self.sp[k] = (CubicSpline(eps, v.real), CubicSpline(eps, v.imag))
            else:
                self.sp[k] = CubicSpline(eps, v)
        self.eps_max = eps_max

    def __call__(self, key, e):
        s = self.sp[key]
        if isinstance(s, tuple):
            return s[0](e) + 1j * s[1](e)
        return s(e)


def f_p(p, e):
    """(F_p(e) - 1)/e^2 = 4 sin^2(e ln p/2) / (e^2 (p - 2 + p^{i e})), regular at e = 0"""
    lp = math.log(p)
    sinc = np.sinc(e * lp / (2 * math.pi))           # sin(e lp/2)/(e lp/2)
    return (lp ** 2) * sinc ** 2 / (p - 2 + np.exp(1j * e * lp))


def f_p_lead(p, e):
    """leading order in 1/p of (F_p - 1)/e^2: 4 sin^2(e ln p/2)/(e^2 p)"""
    lp = math.log(p)
    return (lp ** 2) * np.sinc(e * lp / (2 * math.pi)) ** 2 / p


def d2lnPsi(B, e):
    """d^2/de^2 ln Psi(e,k) = 2 sum_{p|k} sum_m (ln p)^2 cos(m e ln p)/p^m"""
    out = np.zeros_like(e)
    for p in B:
        lp = math.log(p); m = 1
        while p ** m < 1e18:
            out += lp ** 2 * np.cos(m * e * lp) / p ** m; m += 1
    return 2 * out


def dS_diag(p, L, l):
    """Berry's per-prime diagonal drop, -(2/pi^2) sum_m sin^2(pi l m ln p/L)/(m^2 p^m) (closed form)"""
    lp = math.log(p); m = 1; s = 0.0
    while p ** m < 1e18:
        s += math.sin(math.pi * l * m * lp / L) ** 2 / (m * m * p ** m); m += 1
    return -(2 / math.pi ** 2) * s


# quadrature: pieces of width W from 0 to LMAX, Gauss-Legendre with NG nodes each (period-aligned: 40 pieces per period)
W = 0.025; LMAX = 8.0; NG = 8
_x, _w = np.polynomial.legendre.leggauss(NG)
NPC = int(round(LMAX / W))
R_NODES = (np.arange(NPC)[:, None] * W + (W / 2) * (_x[None, :] + 1)).ravel()
R_W = np.tile((W / 2) * _w, NPC)


def cum_sigma(fvals):
    """given f on R_NODES, return S(l) = 2 int_0^l (l - r) f(r) dr at l = W, 2W, ..., LMAX (array of length NPC)"""
    f = fvals.reshape(NPC, NG); wr = R_W.reshape(NPC, NG); r = R_NODES.reshape(NPC, NG)
    I0 = np.cumsum((f * wr).sum(1)); I1 = np.cumsum((f * wr * r).sum(1))
    ls = np.arange(1, NPC + 1) * W
    return 2 * (ls * I0 - I1)


LGRID = np.arange(1, NPC + 1) * W


def sigma2_gue():
    f = -np.sinc(R_NODES) ** 2
    return LGRID + cum_sigma(f)


def C_of_r(G, L, B, r=R_NODES, flip=False):
    """regular part C(r) of R(r)-1 for conductor with bad-prime set B at fixed L. flip=True uses e^{-i eps L}
    (wrong-sign convention, for the sign check)."""
    e = TWO_PI * r / L
    Z2 = G('Z2', e); Z2m1 = G('Z2m1', e); Phi = G('Phi', e); Phim1 = G('Phim1', e)
    Psi = np.ones_like(e, dtype=complex); Psim1 = np.zeros_like(e, dtype=complex)
    for p in B:
        fp = f_p(p, e); F = 1 + e ** 2 * fp
        Psim1 = Psim1 * F + fp; Psi = Psi * F
    rego = Z2m1 * Phi * Psi + Phim1 * Psi + Psim1          # (eps^2|zeta|^2 Phi Psi - 1)/eps^2
    ph = np.exp((-1 if flip else 1) * 1j * TWO_PI * r)
    return (2 / L ** 2) * (G('regd', e) + (ph * rego).real) - (1 / L ** 2) * (G('d2lnPhiDiag', e) + d2lnPsi(B, e))


def sigma2_full(G, L, B, flip=False):
    """absolute BK Sigma2(l) on LGRID at fixed L"""
    return sigma2_gue() + cum_sigma(C_of_r(G, L, B, flip=flip))


def dsig_off(G, L, p, rest=(), lead=False, r=R_NODES):
    """off-diagonal change of Sigma2 from adding p to the bad set `rest` (exact), on LGRID:
    2 int (l-r) (2/L^2) Re[e^{2 pi i r} eps^2|zeta|^2 Phi_off Psi_rest (F_p-1)/eps^2] dr"""
    e = TWO_PI * r / L
    base = G('Z2', e) * G('Phi', e)
    for q in rest:
        base = base * (1 + e ** 2 * f_p(q, e))
    fp = f_p_lead(p, e) if lead else f_p(p, e)
    fv = (2 / L ** 2) * (np.exp(1j * TWO_PI * r) * base * fp).real
    return cum_sigma(fv)


def dsig_diag_grid(p, L):
    return np.array([dS_diag(p, L, l) for l in LGRID])


def L_nodes(q, t1, t2, n=10):
    """Gauss-Legendre nodes in t on [t1, t2], weights proportional to the zero density ln(q t/2pi); returns (L_i, w_i)"""
    x, w = np.polynomial.legendre.leggauss(n)
    t = (t1 + t2) / 2 + (t2 - t1) / 2 * x
    Lt = np.log(q * t / TWO_PI)
    ww = w * Lt
    return Lt, ww / ww.sum()
