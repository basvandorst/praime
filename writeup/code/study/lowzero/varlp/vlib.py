"""vlib.py: two- and three-point functions of the zeros of L(s, chi) (ratios conjecture; bad primes deleted) on a long
range of separations. Same formulas as the repo's bklib.py (two points) and tplib.py (three points); the only change is
that the tables of zeta and of the Euler product Phi are read from xgrid.npz (tab.py, eps up to 160) instead of
grid.npz/gtab.npz (eps up to 7.6). check_vlib.py compares every function with bklib/tplib on the common range.
Unfolded separations r (mean spacing 1), raw eps = 2 pi r / L.
"""
import os, sys, math
import numpy as np
from scipy.interpolate import CubicSpline
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, 'lib'))
import bklib  # noqa: E402  (f_p, primes_upto)

TWO_PI = 2 * math.pi
EG = 0.5772156649015329


class XG:
    def __init__(self, path=os.path.join(HERE, 'lib', 'xgrid.npz')):
        z = np.load(path); e = z['eps']; self.emax = float(e[-1])
        self.sp = {}
        for k in ('Z2', 'Z2m1', 'regd', 'g', 'Phi', 'Phim1', 'd2lnPhiDiag'):
            v = z[k]
            self.sp[k] = (CubicSpline(e, v.real), CubicSpline(e, v.imag)) if np.iscomplexobj(v) else CubicSpline(e, v)

    def __call__(self, k, e):
        a = np.abs(e)
        if np.any(a > self.emax): raise ValueError('eps beyond table: %g' % a.max())
        s = self.sp[k]
        return s[0](a) + 1j * s[1](a) if isinstance(s, tuple) else s(a)


X = XG()


def primes_upto(n): return bklib.primes_upto(n)


# ------------------------------------------------------------------ two points (bklib.C_of_r, extended)
def f_p(p, e): return bklib.f_p(p, np.abs(e))


def T2_gue(r): return -np.sinc(r) ** 2


def C_of_r(L, S, r):
    """T2 - T2_GUE for bad set S (bklib.C_of_r with the extended tables), r >= 0 or any sign (even)."""
    r = np.abs(r); e = TWO_PI * r / L
    Z2m1 = X('Z2m1', e); Phi = X('Phi', e); Phim1 = X('Phim1', e)
    Psi = np.ones_like(e, dtype=complex); Psim1 = np.zeros_like(e, dtype=complex)
    for p in S:
        fp = f_p(p, e); F = 1 + e ** 2 * fp; Psim1 = Psim1 * F + fp; Psi = Psi * F
    rego = Z2m1 * Phi * Psi + Phim1 * Psi + Psim1
    ph = np.exp(1j * TWO_PI * r)
    return (2 / L ** 2) * (X('regd', e) + (ph * rego).real) - (1 / L ** 2) * (X('d2lnPhiDiag', e) + bklib.d2lnPsi(S, e))


def dT2_off(L, p, r, rest=()):
    """off-diagonal change of T2 when p is added to the bad set `rest` (bklib.dsig_off integrand)."""
    r = np.abs(r); e = TWO_PI * r / L
    base = X('Z2', e) * X('Phi', e)
    for q in rest: base = base * (1 + e ** 2 * f_p(q, e))
    return (2 / L ** 2) * (np.exp(1j * TWO_PI * r) * base * f_p(p, e)).real


def dT2_diag(L, p, r):
    e = TWO_PI * np.abs(r) / L; lp = math.log(p); out = np.zeros_like(e); m = 1
    while p ** m < 1e18:
        out += lp ** 2 * np.cos(m * e * lp) / p ** m; m += 1
    return -(2 / L ** 2) * out


# ------------------------------------------------------------------ three points (tplib.T3_parts, extended)
def g_reg(e, S=()):
    v = X('g', e); v = np.where(e < 0, np.conj(v), v)
    for p in S: v = v + math.log(p) / (p ** (1 + 1j * e) - 1)
    return v


def Zc_reg(e, S=()):
    a = np.abs(e)
    Z2m1 = X('Z2m1', a); Phi = X('Phi', a); Phim1 = X('Phim1', a)
    Psi = np.ones_like(a, dtype=complex); Psim1 = np.zeros_like(a, dtype=complex)
    for p in S:
        fp = bklib.f_p(p, a); F = 1 + a ** 2 * fp; Psim1 = Psim1 * F + fp; Psi = Psi * F
    rego = Z2m1 * Phi * Psi + Phim1 * Psi + Psim1
    return np.where(e > 0, np.conj(rego), rego)


class PrimeSums:
    def __init__(self, PB=2000):
        self.P = primes_upto(PB).astype(float); self.lnP = np.log(self.P)

    def B1_all(self, e1, e2, chunk=256):
        out = np.empty(len(e1), complex); P = self.P; lp = self.lnP
        for j in range(0, len(e1), chunk):
            E1 = np.exp(-1j * np.outer(e1[j:j + chunk], lp)); E2 = np.exp(-1j * np.outer(e2[j:j + chunk], lp))
            num = (1 - E1) * (1 - E1 - E2 + E2 / P) * lp
            den = (1 - E2 / (P * E1)) * (1 - E2 / P) * (1 - 2 / P + E1 / P) * P ** 2 * E1 / E2
            out[j:j + chunk] = (num / den).sum(1)
        return out


def b1_p(p, e1, e2):
    lp = math.log(p); E1 = np.exp(-1j * e1 * lp); E2 = np.exp(-1j * e2 * lp)
    return (1 - E1) * (1 - E1 - E2 + E2 / p) * lp / ((1 - E2 / (p * E1)) * (1 - E2 / p) * (1 - 2 / p + E1 / p) * p ** 2 * E1 / E2)


TERMS = (('w', 'u', 'v'), ('v', 'w', 'u'), ('u', 'w', 'v'))


def b1_args(u, v, w, L):
    """the six (f1, f2) argument pairs of B1 in the order used by T3P"""
    k = TWO_PI / L; pos = {'u': u * k, 'v': v * k, 'w': w * k}; out = []
    for (A, B1n, B2n) in TERMS:
        al = -pos[A]; e1 = al + pos[B1n]; e2 = al + pos[B2n]
        out += [(e1, e2), (e2, e1)]
    return out


def T3P(u, v, w, L, Slist, B1all):
    """P-part of the connected T3 (arithmetic P-terms minus their GUE limit, /L^3, six terms, real) for each bad set
    in Slist (list of tuples), given B1all = list of the six all-prime B1 arrays (b1_args order). Returns list."""
    k = TWO_PI / L; pos = {'u': u * k, 'v': v * k, 'w': w * k}
    res = []
    for S in Slist:
        totP = 0; ib = 0
        for (A, B1n, B2n) in TERMS:
            al = -pos[A]; be1 = pos[B1n]; be2 = pos[B2n]
            e1 = al + be1; e2 = al + be2; ed = be2 - be1
            for (f1, f2, fd) in ((e1, e2, ed), (e2, e1, -ed)):
                B = B1all[ib]; ib += 1
                for p in S: B = B - b1_p(p, f1, f2)
                ph = np.exp(-1j * L * f1)
                h = B + g_reg(f2, S) - g_reg(fd, S)
                P3G = -1 / (1j * f2) + 1 / (1j * fd)
                totP = totP - ph * (Zc_reg(f1, S) * (P3G + h) + h / f1 ** 2)
        res.append((2 / L ** 3 * totP).real)
    return res


def T3_gue(u, v, w): return 2 * np.sinc(v - u) * np.sinc(w - v) * np.sinc(w - u)
