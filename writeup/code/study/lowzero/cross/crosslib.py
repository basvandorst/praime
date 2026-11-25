"""Cross-correlation of the zeros of two Dirichlet L-functions: the diagonal (shared-orbit) law and helpers.

Law (derivation.md sections 2-4; ratios recipe = semiclassical diagonal sum, conjecture):
  zeros 1/2 + i gamma of L(s, chi_1) and L(s, chi_2), psi = chi_1 conj(chi_2) (zero at p | k1 k2),
  R12(t1, t2) = dbar_1 dbar_2 + (1/2 pi^2) Re D_psi(-i (t2 - t1)) + [swapped term, section 3],
  D_psi(z) = sum_n Lambda(n)^2 psi(n) n^(-1-z) = (L'/L)'(1+z, psi) - B_psi(z).
  Window counts in the common variable t: dn_j(t) = #{gamma_j in (t, t+w]} - (theta_j(t+w) - theta_j(t))/pi;
  Cov(dn_1, dn_2)(w) = (2/pi^2) sum_p sum_m Re[psi(p)^m] sin^2(m w ln p / 2) / (m^2 p^m)
                     = (1/pi^2) { ln|L(1,psi)| - (ln|L(1+iw,psi)| + ln|L(1-iw,psi)|)/2 } - (1/pi^2) corr2(w).
"""
import math
import numpy as np
from scipy.special import loggamma

PI2 = math.pi ** 2


def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]


def kron_table(Ds, P):
    """chi_D(p) for fundamental discriminants Ds (rows) and primes P (columns), int8 in {-1, 0, 1}."""
    Ds = np.asarray(Ds, np.int64); P = np.asarray(P, np.int64)
    out = np.zeros((len(Ds), len(P)), np.int8)
    for j, p in enumerate(P):
        if p == 2:
            r = Ds % 8
            out[:, j] = np.where(Ds % 2 == 0, 0, np.where((r == 1) | (r == 7), 1, -1))
            continue
        a = Ds % p; e = (p - 1) // 2; res = np.ones_like(a); b = a.copy()
        while e:
            if e & 1: res = (res * b) % p
            b = (b * b) % p; e >>= 1
        out[:, j] = np.where(a == 0, 0, np.where(res == 1, 1, -1))
    return out


def AB(p, w, mmax_val=1e-14):
    """per-prime coefficients of the diagonal covariance for real characters at window widths w (array):
    A_p(w) = (2/pi^2) sum_{m odd} sin^2(m w ln p/2)/(m^2 p^m)  (multiplies chi1(p)chi2(p)),
    B_p(w) = (2/pi^2) sum_{m even} sin^2(m w ln p/2)/(m^2 p^m)  (multiplies 1_{p not | k1 k2})."""
    w = np.asarray(w, float); lp = math.log(p); A = np.zeros_like(w); B = np.zeros_like(w); m = 1
    while True:
        c = 1.0 / (m * m * float(p) ** m)
        if c < mmax_val: break
        t = c * np.sin(m * w * lp / 2) ** 2
        if m % 2: A += t
        else: B += t
        m += 1
    return 2 / PI2 * A, 2 / PI2 * B


def AB_complex(p, w, phase, mmax_val=1e-14):
    """general psi(p) = e^{i phase}: (2/pi^2) sum_m cos(m phase) sin^2(m w ln p/2)/(m^2 p^m)"""
    w = np.asarray(w, float); lp = math.log(p); S = np.zeros_like(w); m = 1
    while True:
        c = 1.0 / (m * m * float(p) ** m)
        if c < mmax_val: break
        S += c * math.cos(m * phase) * np.sin(m * w * lp / 2) ** 2; m += 1
    return 2 / PI2 * S


def theta(t, q, a):
    """theta_chi(t) = (t/2) ln(q/pi) + Im lnGamma((1/2 + a + it)/2); N(t) = theta/pi + const + S(t)"""
    return np.imag(loggamma(0.25 + a / 2 + 0.5j * np.asarray(t, float))) + np.asarray(t, float) / 2 * math.log(q / math.pi)


def dn_series(g, q, a, w, tgrid):
    """window-count fluctuation dn(t) = #{gamma in (t, t+w]} - (theta(t+w) - theta(t))/pi on tgrid"""
    c = np.searchsorted(g, tgrid + w, 'right') - np.searchsorted(g, tgrid, 'right')
    return c - (theta(tgrid + w, q, a) - theta(tgrid, q, a)) / math.pi


def pair_cov_full(X, W, P, Y=None, chunk=4000):
    """predicted diagonal Cov for all pairs of rows of X (chi values at primes P) at widths W:
    returns array (len(W), N, M) with M = rows of Y (default Y = X). Real characters."""
    Y = X if Y is None else Y
    out = np.zeros((len(W), X.shape[0], Y.shape[0]))
    for s in range(0, len(P), chunk):
        Pc = P[s:s + chunk]
        Ac = np.zeros((len(W), len(Pc))); Bc = np.zeros((len(W), len(Pc)))
        for j, p in enumerate(Pc):
            a, b = AB(int(p), W); Ac[:, j] = a; Bc[:, j] = b
        Xc = X[:, s:s + chunk].astype(np.float64); Yc = Y[:, s:s + chunk].astype(np.float64)
        X2 = Xc * Xc; Y2 = Yc * Yc
        for k in range(len(W)):
            out[k] += (Xc * Ac[k]) @ Yc.T + (X2 * Bc[k]) @ Y2.T
    return out
