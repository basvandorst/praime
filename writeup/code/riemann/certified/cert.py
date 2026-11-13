"""Certified explicit-formula samples with a compactly supported kernel.

For a sample point x > 1 and kernel K (support |v| < a), the explicit formula gives exactly

  sum_m Lambda(m) sqrt(x/m) K(x ln(m/x)) = P(x) - (2/sqrt x) sum_{rho} Re(x^{i gamma} Khat(gamma/x)) - E_triv(x)

(sum over all nontrivial zeros; for gamma <= 3e12 they are on the line, Platt-Trudgian 2021).
sample(x, ...) returns the zero side computed from the listed zeros up to T, and a rigorous bound
delta(x) on |zero side computed - prime side|, under

  Hypothesis Z(T, eps): the list holds every zero with 0 < gamma <= T, each within eps.

delta = tail (zeros above T, via Rosser's bound on N(t)) + precision (eps) + trivial zeros + rounding.
"""
import math
import numpy as np
from kbm import Khat, majorant, P, K0, K

def F(t):
    return t / (2 * math.pi) * math.log(t / (2 * math.pi)) - t / (2 * math.pi) + 7 / 8

def Nplus(t):
    # Rosser (1941), Thm 19: |N(T) - F(T)| < 0.137 log T + 0.443 log log T + 1.588 for T >= 2
    return F(t) + 0.137 * math.log(t) + 0.443 * math.log(math.log(t)) + 1.588

def tail_bound(x, T, NT, a, beta, m, steps=4000):
    """Rigorous upper bound of sum_{gamma > T} |Khat(gamma/x)| given exactly NT zeros in (0, T]."""
    assert a * T / x > beta * 1.000001
    p = m + 0.5 + 1 / 3  # majorant decays like (a w)^(-p) for large w
    tK = T * 1e6
    ts = T * np.exp(np.linspace(0, math.log(tK / T), steps + 1))
    f = majorant(ts / x, a, beta, m)
    f = np.minimum.accumulate(f)  # enforce non-increasing (it already is)
    Np = ts / (2 * math.pi) * np.log(ts / (2 * math.pi)) - ts / (2 * math.pi) + 7 / 8 + 0.137 * np.log(ts) + 0.443 * np.log(np.log(ts)) + 1.588
    s = float(np.sum((Np[1:] - NT) * (f[:-1] - f[1:])))
    # remainder beyond tK: int_{tK}^inf (N+(t) - NT)(-f') dt <= N+(tK) f(tK) + int_{tK}^inf f(t) dN+(t)
    # f(t) <= C t^-p for t >= tK with C = f(tK) tK^p; dN+/dt <= log(t)/(2 pi) + 1/t
    C = f[-1] * tK ** p
    integral = C * (tK ** (1 - p)) * ((math.log(tK) / (2 * math.pi)) / (p - 1) + 1 / (2 * math.pi * (p - 1) ** 2) + 1 / p)
    return s + Nplus(tK) * f[-1] + integral

H_RH = 3e12   # RH verified for 0 < gamma <= 3e12 (Platt and Trudgian, Bull. LMS 53 (2021))
_C2 = {}
def offline_bound(x, a, beta, m):
    """Bound on the whole contribution of zeros with gamma > H_RH, on or off the critical line.

    sqrt(x)|phi~(rho)| <= C2 x^2 / gamma^2 by two integrations by parts (K, K' vanish at +-a for m >= 2),
    C2 = e^{a/(2x)} int (|K''| + |K'|/x + |K|/(4x^2)) dv; summed over gamma > H with
    sum_{gamma > H} gamma^-2 <= int_H^inf 2 N+(t) t^-3 dt <= 1.1 ln(H) / (pi H).
    """
    assert m >= 2
    if m >= 3:
        return offline_bound_j(x, a, beta, m, min(m, 4))
    key = (a, beta, m)
    if key not in _C2:
        v = np.linspace(-a, a, 200001); h = v[1] - v[0]
        k = K(v, a, beta, m); k1 = np.gradient(k, h); k2 = np.gradient(k1, h)
        _C2[key] = (float(np.sum(np.abs(k2))) * h * 1.05, float(np.sum(np.abs(k1))) * h * 1.05, float(np.sum(np.abs(k))) * h * 1.05)
    i2, i1, i0 = _C2[key]
    C2 = math.exp(a / (2 * x)) * (i2 + i1 / x + i0 / (4 * x * x))
    return 2 * C2 * x * x * 1.1 * math.log(H_RH) / (math.pi * H_RH)

_CJ = {}
def offline_bound_j(x, a, beta, m, j):
    """Same bound with j integrations by parts in v (allowed for j <= m: K^(i)(+-a) = 0 for i < m).

    phi~(rho) = x^(rho - 1) int K(v) e^{s v/x} dv with s = rho - 1/2, |Re s| < 1/2, |s| >= gamma, so
    |phi~(rho)| <= x^(sigma - 1) (x/gamma)^j e^{a/(2x)} int |K^(j)| <= e^{a/(2x)} I_j x^j gamma^-j.
    Zeros with gamma > H_RH, each with its conjugate: <= 2 e^{a/(2x)} I_j x^j sum_{gamma > H} gamma^-j, and
    sum_{gamma > H} gamma^-j <= int_H^inf N+(t) j t^(-j-1) dt <= (j / 2 pi) H^(1-j) (ln H/(j-1) + 1/(j-1)^2)
    with N+(t) <= t ln(t) / (2 pi) for t >= H.
    """
    key = (a, beta, m, j)
    if key not in _CJ:
        v = np.linspace(-a, a, 400001); h = v[1] - v[0]
        k = K(v, a, beta, m)
        for _ in range(j):
            k = np.gradient(k, h)
        _CJ[key] = float(np.sum(np.abs(k))) * h * 1.1
    I = _CJ[key]
    H = H_RH
    s = (j / (2 * math.pi)) * H ** (1 - j) * (math.log(H) / (j - 1) + 1 / (j - 1) ** 2)
    return 2 * math.exp(a / (2 * x)) * I * x ** j * s

def sample(x, G, T, a, beta, m, eps):
    """Zero side at x from zeros G[G <= T], and the rigorous error bound delta."""
    assert T <= G[-1], 'T beyond the last listed zero: the list is not known to be complete there'
    g = G[G <= T]
    NT = len(g)
    gf = g.astype(float)
    kh = Khat(gf / x, a, beta, m)
    L = math.log(x)
    if G.dtype == np.longdouble:   # zeros and phases in 64-bit-mantissa arithmetic
        cg = np.cos(g * np.log(np.longdouble(x))).astype(float); ua = 4.5e-19
    else:
        cg = np.cos(g * L); ua = 4.5e-16
    z = float(P(x, a, beta, m) - (2 / math.sqrt(x)) * np.sum(cg * kh))
    akh = np.abs(kh)
    tail = (2 / math.sqrt(x)) * tail_bound(x, T, NT, a, beta, m)
    # precision: |d/dgamma [cos(gamma L) Khat(gamma/x)]| <= L |Khat| + |Khat'|/x, |Khat'| <= a * Khat(0)
    kh0 = float(Khat(0.0, a, beta, m))
    prec = (2 / math.sqrt(x)) * eps * (L * float(np.sum(akh)) * 1.01 + NT * a * kh0 / x + NT * eps * L * L * kh0)
    # trivial zeros: sqrt(x) * int phi(t)/(t^3 - t) dt with t = x e^{v/x}, computed by quadrature
    assert x * math.exp(-a / x) > 1.0, 'kernel support must lie in t > 1'
    from scipy.integrate import quad
    def fv(v):
        t = x * math.exp(v / x)
        return float(K(v, a, beta, m)) * math.sqrt(t) / (x * (t ** 3 - t))
    tv, terr = quad(fv, -a, a, limit=200, epsabs=0, epsrel=1e-12)
    z -= math.sqrt(x) * tv
    triv = 10 * math.sqrt(x) * terr + 1e-15 * abs(math.sqrt(x) * tv)
    # rounding: argument gamma*L has absolute error <= 4 ulp(gamma L); each term relative error <= 1e-15
    rnd = (2 / math.sqrt(x)) * float(np.sum(akh * (gf * L * ua + 1e-15))) + 1e-14 * P(x, a, beta, m)
    off = offline_bound(x, a, beta, m)
    return z, tail + prec + triv + rnd + off, dict(NT=NT, tail=tail, prec=prec, triv=triv, rnd=rnd, off=off)
