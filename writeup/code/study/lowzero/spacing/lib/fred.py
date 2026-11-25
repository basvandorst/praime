"""Fredholm-determinant tools (Bornemann's Nystrom method) for translation-invariant kernels on an interval.
E(0; s) = det(I - K_[0,s]); spacing density p(s) = E''(s); Var(s) = 2 int_0^inf E(0;s) ds - 1 (unit density).
Janossy density for two consecutive spacings: j(s1, s2) = det(I - K_J) det[R_J(x_a, x_b)]_{a,b in {0, s1, s1+s2}},
J = [0, s1+s2], R_J = K (I - K_J)^-1 (resolvent kernel), by Nystrom interpolation."""
import numpy as np


def gl(m, a, b):
    x, w = np.polynomial.legendre.leggauss(m)
    return (a + b) / 2 + (b - a) / 2 * x, (b - a) / 2 * w


def sinc_kernel(r):
    return np.sinc(r)


def E0(Kfun, s, m=40):
    """det(I - K) on [0, s]; Kfun(r) even function of the separation."""
    if s <= 0: return 1.0
    x, w = gl(m, 0.0, s); sw = np.sqrt(w)
    A = np.eye(m) - sw[:, None] * Kfun(x[:, None] - x[None, :]) * sw[None, :]
    return float(np.linalg.det(A))


def E0_grid(Kfun, smax=6.0, ds=0.01, m=40):
    s = np.arange(0, smax + ds / 2, ds)
    return s, np.array([E0(Kfun, v, m) for v in s])


def var_from_E(s, E):
    """Var = 2 int E ds - 1 (trapezoid on a fine grid; E decays like exp(-pi^2 s^2/8))"""
    return 2 * np.trapezoid(E, s) - 1


def janossy3(Kfun, s1, s2, m=40):
    """density of: points at 0, s1, s1+s2 and none else in [0, s1+s2] (unit-density stationary determinantal)."""
    T = s1 + s2
    x, w = gl(m, 0.0, T); sw = np.sqrt(w)
    Kxx = Kfun(x[:, None] - x[None, :])
    A = np.eye(m) - sw[:, None] * Kxx * sw[None, :]
    det = np.linalg.det(A)
    pts = np.array([0.0, s1, T])
    Kpx = Kfun(pts[:, None] - x[None, :])                  # 3 x m
    Kpp = Kfun(pts[:, None] - pts[None, :])
    # R(a,b) = K(a,b) + sum_k K(a,x_k) w_k f_k(b), with (I - K W) f(.,b) = K(x, b)
    rhs = (Kpx.T) * 1.0                                     # m x 3 : K(x_k, b)
    # solve (I - W^1/2 K W^1/2) g = W^1/2 K(x,b); then f = W^-1/2 g  -> sum_k K(a,x_k) w_k f_k = sum_k K(a,x_k) sw_k g_k
    g = np.linalg.solve(A, sw[:, None] * rhs)
    R = Kpp + (Kpx * sw[None, :]) @ g
    return det * np.linalg.det(R)


def janossy_table(Kfun, nT=40, nf=24, Tmax=5.0, m=28):
    """Janossy density on the (T, f) grid, T = s1 + s2 in (0, Tmax), f = s1/T in (0, 1/2) (Gauss in each);
    by symmetry the region f > 1/2 is the mirror image.  Returns T, f, weights (incl. Jacobian T and the factor 2), j."""
    xT, wT = gl(nT, 0.0, Tmax); xf, wf = gl(nf, 0.0, 0.5)
    J = np.empty((nT, nf))
    for a, T in enumerate(xT):
        for b, f in enumerate(xf):
            J[a, b] = janossy3(Kfun, f * T, (1 - f) * T, m)
    W = 2 * (wT * xT)[:, None] * wf[None, :]
    return xT, xf, W, J


def ratio_mean(Kfun, **kw):
    """<r~> = <min/max> of the two spacings around a typical point, and the normalisation of the Janossy density"""
    xT, xf, W, J = janossy_table(Kfun, **kw)
    r = (xf / (1 - xf))[None, :]
    norm = (W * J).sum()
    return float((W * J * r).sum() / norm), float(norm)
