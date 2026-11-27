"""Measure F_t(x) = sum_k w(gamma_k(t)) exp(i gamma_k(t) ln x) and compare with the frozen prediction
   P_t(n) = -(1/2pi) ln(n) kappa_t(n) n^(-1/2) int w(g) exp(-(t/4) ln(g/2pi) ln n) dg  (times exp(i s ln n)).
w = Kaiser(beta) taper on fixed heights [g_lo, g_hi] (window fraction of the moving zeros)."""
import numpy as np, sys
from scipy.special import i0

BETA = 10.0


def dirlog(t, N):
    """d(n) = ln(n) kappa_t(n): coefficients of -(d/ds) log sum_m exp((t/4) ln^2 m) m^-s."""
    ln = np.log(np.arange(N + 1.0)); ln[0] = 0
    a = np.exp(t / 4 * ln ** 2)
    d = np.zeros(N + 1)
    for n in range(2, N + 1):
        d[n] += a[n] * ln[n]
        for m in range(2 * n, N + 1, n):         # n divides m, n < m: subtract a(m/n) d(n)
            d[m] -= a[m // n] * d[n]
    return d


def kaiser(g, lo, hi, beta=BETA):
    u = (2 * g - lo - hi) / (hi - lo)
    w = np.zeros_like(g); ins = np.abs(u) < 1
    w[ins] = i0(beta * np.sqrt(1 - u[ins] ** 2)) / i0(beta)
    return w


def Fsum(g, w, xs, chunk=200):
    out = np.empty(len(xs), complex)
    for a in range(0, len(xs), chunk):
        om = np.log(xs[a:a + chunk])
        out[a:a + chunk] = (w[:, None] * np.exp(1j * g[:, None] * om[None, :])).sum(0)
    return out


def template(t, ns, lo, hi, dn):
    gg = np.linspace(lo, hi, 20001); ww = kaiser(gg, lo, hi)
    L = np.log(gg / (2 * np.pi)); ln = np.log(ns)
    I = np.trapezoid(ww[None, :] * np.exp(-(t / 4) * L[None, :] * ln[:, None]), gg, axis=1)
    return -(1 / (2 * np.pi)) * dn[ns] * ns ** -0.5 * I


def window(x0, f0, f1):
    M = len(x0); return x0[int(f0 * M)], x0[int(f1 * M) - 1]


def classify(n):
    f = []; m = n; p = 2
    while p * p <= m:
        while m % p == 0: f.append(p); m //= p
        p += 1
    if m > 1: f.append(m)
    if len(f) == 1: return 'prime'
    if len(set(f)) == 1: return 'ppow'
    return 'comp'
