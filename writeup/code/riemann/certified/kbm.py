"""Generalized Kaiser-Bessel kernel of order m (Lewitt's "blob"), compact support |v| < a.

K(v)   = (1 - (v/a)^2)^(m/2) * I_m(beta*sqrt(1 - (v/a)^2))      for |v| < a, else 0
Khat(w)= int K(v) e^{i w v} dv
       = a*sqrt(2 pi)*beta^m * Z(a w),  nu = m + 1/2,
  Z(s) = I_nu(sqrt(beta^2 - s^2)) / (beta^2 - s^2)^(nu/2)   for |s| < beta
       = J_nu(sqrt(s^2 - beta^2)) / (s^2 - beta^2)^(nu/2)   for |s| > beta
K vanishes like (a - |v|)^m at the edges, so Khat decays like |w|^(-m-1).
Rigorous majorant for |s| > beta: |J_nu(y)| <= min(y^nu / (2^nu Gamma(nu+1)), 0.7858 y^(-1/3))  (Landau 2000).
"""
import math
import numpy as np
from scipy.special import iv, jv, ive, gammaln

def K(v, a, beta, m):
    v = np.asarray(v, dtype=float)
    t = np.clip(1 - (v / a) ** 2, 0, None)
    out = t ** (m / 2) * iv(m, beta * np.sqrt(t))
    return np.where(np.abs(v) < a, out, 0.0)

def K0(beta, m):
    return float(iv(m, beta))

def Z(s, beta, m):
    nu = m + 0.5
    s = np.abs(np.asarray(s, dtype=float))
    out = np.empty_like(s)
    lo = s < beta
    q = np.sqrt(beta * beta - s[lo] ** 2)
    small = q < 1e-8
    zl = np.where(small, 1.0, iv(nu, np.where(small, 1.0, q)) / np.where(small, 1.0, q) ** nu)
    zl = np.where(small, 1.0 / (2 ** nu * math.gamma(nu + 1)), zl)
    out[lo] = zl
    y = np.sqrt(s[~lo] ** 2 - beta * beta)
    smallh = y < 1e-8
    zh = np.where(smallh, 1.0, jv(nu, np.where(smallh, 1.0, y)) / np.where(smallh, 1.0, y) ** nu)
    zh = np.where(smallh, 1.0 / (2 ** nu * math.gamma(nu + 1)), zh)
    out[~lo] = zh
    return out

def Khat(w, a, beta, m):
    return a * math.sqrt(2 * math.pi) * beta ** m * Z(a * np.asarray(w, dtype=float), beta, m)

def Khat_deriv_bound(w, a, beta, m, h=1e-5):
    """|d Khat / dw| by a central difference (used only inside the precision term, with a safety factor)."""
    return np.abs(Khat(np.asarray(w) + h, a, beta, m) - Khat(np.asarray(w) - h, a, beta, m)) / (2 * h)

def P(x, a, beta, m):
    """int K(v) e^{v/(2x)} dv = Khat(-i/(2x))."""
    nu = m + 0.5
    r = a / (2 * x)
    q = math.sqrt(beta * beta + r * r)
    return a * math.sqrt(2 * math.pi) * beta ** m * float(iv(nu, q)) / q ** nu

def majorant(w, a, beta, m):
    """Non-increasing M(w) >= |Khat(w')| for all w' >= w, valid for a*w > beta."""
    nu = m + 0.5
    s = a * np.asarray(w, dtype=float)
    y = np.sqrt(np.clip(s * s - beta * beta, 1e-300, None))
    c0 = 1.0 / (2 ** nu * math.gamma(nu + 1))
    b = np.minimum(c0, 0.7858 * y ** (-1.0 / 3.0 - nu))
    return a * math.sqrt(2 * math.pi) * beta ** m * b
