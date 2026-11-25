"""Gap statistics of a re-mapped (deformed) process, from the base's two-point Janossy densities.
Base: the sine process (sinc kernel), tabulated once (jtab.npz).
  p(u)      = density of the spacing following a typical point
  j(u1, u2) = joint density of the spacings before and after a typical point
Deformation by one prime's wave (closure.Wave), phase phi uniform, x0 = 0 (the phase absorbs the position):
  <F>_with = E_phi[ N'(0) * int int F(s1, s2) j(u1(s1), u2(s2)) N'(-s1) N'(s2) ds1 ds2 ],
  u1 = N(0) - N(-s1),  u2 = N(s2) - N(0)   (Palm weighting N'(0); Jacobians N'(-s1), N'(s2)).
Second-order formulas for a stationary deformation with increment structure function v(u) (any sign, also formal):
  dVar  = int v(u) p(u) du
  d<r~> = int int_{u1<u2} [ -(v(u1+u2) - v(u1) - v(u2))/u2^2 + 2 u1 v(u2)/u2^3 ] j du1 du2 - int (4 v(u) - v(2u)) j(u,u)/u du
"""
import os, math
import numpy as np
from scipy.interpolate import RectBivariateSpline, CubicSpline
import fred

HERE = os.path.dirname(os.path.abspath(__file__))


def build_tables(h=0.04, umax=6.0, m=32):
    f = HERE + '/jtab.npz'
    if os.path.exists(f):
        z = np.load(f); return z['u'], z['J'], z['p']
    u = np.arange(0, umax + h / 2, h); n = len(u); J = np.zeros((n, n))
    for a in range(n):
        for b in range(a, n):
            if u[a] + u[b] > 6.4 or u[a] == 0 or u[b] == 0: continue
            J[a, b] = J[b, a] = fred.janossy3(fred.sinc_kernel, u[a], u[b], m)
    # one-gap density p(u) from E'' : Janossy with two points: det(I-K) det[R(0,0) R(0,u); R(u,0) R(u,u)]
    p = np.zeros(n)
    for a in range(1, n):
        x, w = fred.gl(m, 0.0, u[a]); sw = np.sqrt(w)
        A = np.eye(m) - sw[:, None] * np.sinc(x[:, None] - x[None, :]) * sw[None, :]
        pts = np.array([0.0, u[a]]); Kpx = np.sinc(pts[:, None] - x[None, :]); Kpp = np.sinc(pts[:, None] - pts[None, :])
        g = np.linalg.solve(A, sw[:, None] * Kpx.T); R = Kpp + (Kpx * sw[None, :]) @ g
        p[a] = np.linalg.det(A) * np.linalg.det(R)
    np.savez(f, u=u, J=J, p=p)
    return u, J, p


class Base:
    def __init__(self):
        self.u, self.J, self.p = build_tables()
        self.js = RectBivariateSpline(self.u, self.u, self.J, kx=3, ky=3)
        self.ps = CubicSpline(self.u, self.p)
        self.umax = self.u[-1]

    def j(self, u1, u2):
        out = self.js.ev(np.clip(u1, 0, self.umax), np.clip(u2, 0, self.umax))
        out[(u1 < 0) | (u2 < 0) | (u1 + u2 > 6.4)] = 0.0
        return np.maximum(out, 0.0)

    def pdf(self, u):
        out = self.ps(np.clip(u, 0, self.umax)); out[(u < 0) | (u > self.umax)] = 0.0
        return np.maximum(out, 0.0)


def grids(nT=48, nf=20, Tmax=6.2):
    xT, wT = fred.gl(nT, 0.0, Tmax)
    xa, wa = fred.gl(nf, 0.0, 0.5); xb, wb = fred.gl(nf, 0.5, 1.0)
    xf = np.concatenate([xa, xb]); wf = np.concatenate([wa, wb])
    T = np.repeat(xT, len(xf)); f = np.tile(xf, nT); W = np.repeat(wT * xT, len(xf)) * np.tile(wf, nT)
    return T * f, T * (1 - f), W


def stats_with(base, wave, nphi=64, gr=None):
    """(<s>, Var(s), <r~>, normalisation) of the re-mapped process; wave=None gives the base itself"""
    s1, s2, W = grids() if gr is None else gr
    s1d = np.linspace(1e-6, 6.0, 1201)
    rt = np.minimum(s1, s2) / np.maximum(s1, s2)
    if wave is None:
        jj = base.j(s1, s2); nrm = (W * jj).sum()
        pp = base.pdf(s1d); m1 = np.trapezoid(s1d * pp, s1d); m2 = np.trapezoid(s1d ** 2 * pp, s1d)
        return m1, m2 - m1 ** 2, (W * jj * rt).sum() / nrm, nrm
    acc = np.zeros(5)
    for phi in 2 * math.pi * np.arange(nphi) / nphi:
        r0 = wave.rho(0.0, phi); N0 = wave.N(0.0, phi)
        u1 = N0 - wave.N(-s1, phi); u2 = wave.N(s2, phi) - N0
        jj = base.j(u1, u2) * wave.rho(-s1, phi) * wave.rho(s2, phi) * r0
        u = wave.N(s1d, phi) - N0; pp = base.pdf(u) * wave.rho(s1d, phi) * r0
        acc += [(W * jj).sum(), (W * jj * rt).sum(), np.trapezoid(pp, s1d), np.trapezoid(s1d * pp, s1d), np.trapezoid(s1d ** 2 * pp, s1d)]
    acc /= nphi
    m1 = acc[3] / acc[2]; m2 = acc[4] / acc[2]
    return m1, m2 - m1 ** 2, acc[1] / acc[0], acc[0]


def second_order(base, v):
    """dVar and d<r~> for a deformation with structure function v (callable, even, v(0) = 0), base = sine."""
    uu = np.linspace(1e-6, 6.0, 3001)
    dvar = np.trapezoid(v(uu) * base.pdf(uu), uu)
    # region u1 < u2: u1 = f T, u2 = (1-f) T, f in (0, 1/2)
    xT, wT = fred.gl(60, 0.0, 6.3); xf, wf = fred.gl(30, 0.0, 0.5)
    T = np.repeat(xT, len(xf)); f = np.tile(xf, len(xT)); W = np.repeat(wT * xT, len(xf)) * np.tile(wf, len(xT))
    u1, u2 = f * T, (1 - f) * T
    g = -(v(u1 + u2) - v(u1) - v(u2)) / u2 ** 2 + 2 * u1 * v(u2) / u2 ** 3
    reg = (W * g * base.j(u1, u2)).sum()
    ud = np.linspace(1e-4, 3.15, 1500)
    kink = np.trapezoid((4 * v(ud) - v(2 * ud)) * base.j(ud, ud) / ud, ud)
    return dvar, reg - kink
