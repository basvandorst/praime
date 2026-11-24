"""Independent implementation of Sigma^2(l) from the Bogomolny-Keating R2 (arXiv:1307.6010).

Sigma^2(l) = Sigma^2_GUE(l)
           - (1/2pi^2) [G(dE) - G(0)],  G = ln|eps zeta(1+i eps)|^2 + ln Phi_diag + ln Psi_diag(k)   (diag, exact double integral)
           + 2 int_0^l (l-r) (1/(2 pi^2 r^2)) * 2? ...  see off_int below
with dE = 2 pi l / L, eps = 2 pi r / L, L = ln(k t / 2 pi).
Off-diagonal: R2off/dbar^2 = (2/L^2) Re[e^{2 pi i r} |zeta|^2 Phi_off Psi_off]
  = (1/(2 pi^2 r^2)) Re[e^{2 pi i r} X(eps)],  X = |eps zeta|^2 Phi_off Psi_off   (since 2/L^2 * 1/eps^2 = 1/(2 pi^2 r^2))
GUE part uses X -> 1 together with the -2 ln|eps| piece of the diagonal.
"""
import numpy as np, mpmath
from scipy.interpolate import CubicSpline
from sympy import primerange, factorint

PMAX = 200000
PR = np.array(list(primerange(2, PMAX)), dtype=float)
LNP = np.log(PR)

EPS = np.linspace(0.0, 8.0, 8001)  # eps grid step 1e-3

def build_tables():
    # ln|eps zeta(1+i eps)|^2
    lnh = np.zeros_like(EPS)
    for i, e in enumerate(EPS):
        if e == 0:
            lnh[i] = 0.0
        else:
            z = mpmath.zeta(mpmath.mpc(1, e))
            lnh[i] = float(mpmath.log(abs(e * z) ** 2))
    # ln Phi_diag = 2 sum_p sum_{m>=2} (1-m)/(m^2 p^m) cos(m eps ln p)
    lnPhiD = np.zeros_like(EPS)
    # ln Phi_off = sum_p ln(1 - (p^{i eps}-1)^2/(p-1)^2)
    lnPhiO = np.zeros_like(EPS, dtype=complex)
    for chunk in np.array_split(np.arange(len(PR)), 200):
        p = PR[chunk][:, None]; lp = LNP[chunk][:, None]
        x = np.exp(1j * EPS[None, :] * lp)
        fac = 1 - (x - 1) ** 2 / (p - 1) ** 2
        fac[p[:, 0] < 4] = 1.0   # p = 2, 3 handled exactly (p = 3 factor has a zero at eps = pi/ln 3)
        lnPhiO += np.log(fac).sum(0)
        for m in range(2, 60):
            w = (1 - m) / (m * m * p ** m)
            if np.all(np.abs(w) < 1e-18):
                break
            lnPhiD += (2 * w * np.cos(m * EPS[None, :] * lp)).sum(0)
    return lnh, lnPhiD, lnPhiO

def lnPsiD_p(eps, p):
    s = np.zeros_like(eps)
    for m in range(1, 80):
        if p ** -m < 1e-18:
            break
        s += -2 * np.cos(m * eps * np.log(p)) / (m * m * p ** m)
    return s

def F_p(eps, p):
    x = p ** (1j * eps)
    return (p - 1 / x) / (p - 2 + x)

class Pred:
    def __init__(self, tables=None):
        if tables is None:
            tables = build_tables()
        self.lnh, self.lnPhiD, self.lnPhiO = tables
        self.s_lnh = CubicSpline(EPS, self.lnh)
        self.s_lnPhiD = CubicSpline(EPS, self.lnPhiD)
        self.s_PhiO_re = CubicSpline(EPS, self.lnPhiO.real)
        self.s_PhiO_im = CubicSpline(EPS, self.lnPhiO.imag)
        # GUE on fine r grid
        self.r = np.linspace(0, 8.0, 16001)
        self.dr = self.r[1] - self.r[0]

    def sigma2(self, ls, L, badp, mode='full', offsign=1.0, prodS=True, Lphase=None):
        """Sigma^2 at lengths ls for density log L, bad primes badp.
        mode: 'full' (BK), 'diag' (bad primes only in diagonal), 'zeta' (no bad primes at all)
        """
        ls = np.asarray(ls, float)
        r = self.r
        eps = 2 * np.pi * r / L
        # GUE number variance via integral
        with np.errstate(divide='ignore', invalid='ignore'):
            gue_int = np.where(r > 0, -np.sin(np.pi * r) ** 2 / (np.pi ** 2 * r ** 2), -1.0)
        X = np.exp(self.s_lnh(eps) + self.s_PhiO_re(eps) + 1j * self.s_PhiO_im(eps))
        for p in (2.0, 3.0):
            x = p ** (1j * eps)
            if mode == 'full' and int(p) in badp:
                X = X * (p - x) * (p - 1 / x) / (p - 1) ** 2   # Phi_off factor times F_p, pole cancelled
            else:
                X = X * (1 - (x - 1) ** 2 / (p - 1) ** 2)
        if mode == 'full':
            for p in badp:
                if p > 3:
                    X = X * F_p(eps, p)
        phase = 2 * np.pi * r if Lphase is None else eps * Lphase
        with np.errstate(divide='ignore', invalid='ignore'):
            off_int = np.where(r > 0, offsign * np.real(np.exp(1j * phase) * (X - 1)) / (2 * np.pi ** 2 * r ** 2), 0.0)
        # fix r=0 by extrapolation
        off_int[0] = 2 * off_int[1] - off_int[2]
        g = gue_int + off_int  # connected R2 (unfolded), beyond GUE's 2-term
        # cumulative integrals int_0^l g and int_0^l r g (trapezoid, fine grid)
        def cumtrap(y):
            c = np.concatenate([[0], np.cumsum((y[1:] + y[:-1]) * 0.5 * self.dr)])
            return c
        G1 = cumtrap(g); G2 = cumtrap(r * g)
        G1l = np.interp(ls, r, G1); G2l = np.interp(ls, r, G2)
        s2 = ls + 2 * (ls * G1l - G2l)
        dE = 2 * np.pi * ls / L
        Gd = self.s_lnh(dE) + self.s_lnPhiD(dE) - (0 + float(self.s_lnPhiD(0.0)))
        if mode in ('full', 'diag'):
            for p in badp:
                Gd = Gd + lnPsiD_p(dE, p) - lnPsiD_p(np.zeros(1), p)
        s2 = s2 - Gd / (2 * np.pi ** 2)
        return s2

def badprimes(k):
    return sorted(factorint(abs(k)).keys())

def band_average(pred, ls, k, t0, t1, nodes=16, **kw):
    x, w = np.polynomial.legendre.leggauss(nodes)
    t = 0.5 * (t1 - t0) * x + 0.5 * (t1 + t0)
    w = w * 0.5 * (t1 - t0)
    Ls = np.log(abs(k) * t / (2 * np.pi))
    wt = w * Ls
    bp = badprimes(k)
    acc = 0
    for Li, wi in zip(Ls, wt):
        acc = acc + wi * pred.sigma2(ls, Li, bp, **kw)
    return acc / wt.sum()
