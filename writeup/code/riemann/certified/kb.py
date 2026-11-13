"""Compactly supported Kaiser-Bessel kernel for the explicit formula.

K(v) = I0(beta*sqrt(1 - (v/a)^2)) - 1 on |v| < a, 0 outside.  K is continuous, vanishes at +-a.
Khat(w) = int K(v) e^{i w v} dv  (real, even), closed form.
Identity (explicit formula, phi(t) = t^{-1/2} K(x ln(t/x)) / ... see paper):

  sum_m Lambda(m) sqrt(x/m) K(x ln(m/x)) = P(x) - (2/sqrt x) sum_{gamma>0} Re(x^{i gamma} Khat(gamma/x)) - triv(x)

with P(x) = int K(v) e^{v/(2x)} dv.
"""
import numpy as np

def K(v, a, beta):
    v = np.asarray(v, dtype=float)
    t = np.clip(1 - (v / a) ** 2, 0, None)
    out = np.i0(beta * np.sqrt(t)) - 1.0
    return np.where(np.abs(v) < a, out, 0.0)

def Khat(w, a, beta):
    """Fourier transform of K at real w (vectorised)."""
    w = np.abs(np.asarray(w, dtype=float)) * a
    out = np.empty_like(w)
    sw = np.where(w > 0, np.sin(w) / np.where(w > 0, w, 1), 1.0)
    lo = w < beta
    q = np.sqrt(np.where(lo, beta * beta - w * w, 0.0))
    out[lo] = np.where(q[lo] > 1e-12, np.sinh(q[lo]) / np.where(q[lo] > 1e-12, q[lo], 1), 1.0) - sw[lo]
    hi = ~lo
    s = np.sqrt(np.where(hi, w * w - beta * beta, 0.0))
    out[hi] = np.where(s[hi] > 1e-12, np.sin(s[hi]) / np.where(s[hi] > 1e-12, s[hi], 1), 1.0) - sw[hi]
    return 2 * a * out

def Khat_prime(w, a, beta, h=1e-6):
    return (Khat(w + h, a, beta) - Khat(w - h, a, beta)) / (2 * h)

def P(x, a, beta):
    """int K(v) e^{v/(2x)} dv in closed form."""
    r = a / (2 * x)
    q = np.sqrt(beta * beta + r * r)
    return 2 * a * (np.sinh(q) / q - np.sinh(r) / r)

def K0(beta):
    return float(np.i0(beta) - 1)
