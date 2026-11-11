"""The explicit formula backwards: zeta zeros from the prime powers up to X.

N_X(t) = theta(t)/pi + 1 - (1/pi) sum_{n <= X} Lambda(n) w(ln n / ln X) sin(t ln n) / (sqrt(n) ln n)
with w = 1 on [0, 6/7] and a cosine fade to 0 at 1. By Weil's explicit formula
N_X(t) = C + sum_gamma Kc(ln X (t - gamma)) + pole term, where
Kc(s) = 1/2 + (1/pi) int_0^1 w(v) sin(s v) / v dv is a smoothed unit step and the pole term
(2/pi) int_0^U w(u/U) cosh(u/2) sin(t u)/u du (U = ln X) tends to 1.
A least-squares fit of the step positions gives the zeros (deconvolution).
Assumes the zeros are simple and on the critical line.

usage:
  python3 primes_to_zeros.py threshold true.f64        windows of 60 zeros at T = 2 pi (X / r)^2
  python3 primes_to_zeros.py all X TMAX KEEP true.f64 [out.f64]
                                                      every zero up to KEEP from the primes up to X
true.f64: zeros from zeros.js, used only to score the result.
"""
import sys
import time

import numpy as np
from scipy.optimize import least_squares
from scipy.special import sici

A = 6 / 7


def wfun(u):
    return np.where(u < A, 1.0, 0.5 * (1 + np.cos(np.pi * (u - A) * 7)))


def lam_table(N):
    s = np.ones(N + 1, bool)
    s[:2] = False
    for p in range(2, int(N ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    ps = np.nonzero(s)[0]
    ns, ls = [ps], [np.log(ps)]
    k = 2
    while True:
        q = ps[ps.astype(float) ** k <= N]
        if len(q) == 0:
            break
        ns.append(q ** k)
        ls.append(np.log(q))
        k += 1
    n = np.concatenate(ns).astype(np.float64)
    o = np.argsort(n)
    return n[o], np.concatenate(ls)[o]


def theta(t):
    return t / 2 * np.log(t / (2 * np.pi)) - t / 2 - np.pi / 8 + 1 / (48 * t) + 7 / (5760 * t ** 3)


H, SM = 0.002, 400.0
S = np.arange(-SM, SM + H / 2, H)
xg, wg = np.polynomial.legendre.leggauss(256)
V = A + (xg + 1) / 2 * (1 - A)
WV = wg * (1 - A) / 2 * 0.5 * (1 + np.cos(np.pi * (V - A) * 7))
K1, KC = np.empty_like(S), np.empty_like(S)
for i in range(0, len(S), 20000):
    s = S[i:i + 20000]
    ss = np.where(np.abs(s) < 1e-12, 1e-12, s)
    K1[i:i + 20000] = (np.where(np.abs(s) < 1e-12, A, np.sin(A * ss) / ss) + np.cos(np.outer(s, V)) @ WV) / np.pi
    KC[i:i + 20000] = 0.5 + (np.sign(s) * sici(A * np.abs(s))[0] + (np.sin(np.outer(s, V)) / V) @ WV) / np.pi


def k1(s):
    return np.interp(s, S, K1, left=0, right=0)


def kc(s):
    return np.interp(s, S, KC, left=0, right=1)


class Primes:
    def __init__(self, X):
        n, l = lam_table(int(X))
        self.U = np.log(X)
        self.ln = np.log(n)
        self.coef = l * wfun(self.ln / self.U) / (np.sqrt(n) * self.ln)

    def pole(self, t):
        out = np.ones_like(t)
        nu = 200000
        u = (np.arange(nu) + 0.5) * self.U / nu
        f = wfun(u / self.U) * np.cosh(u / 2) / u * self.U / nu
        for i in np.nonzero(t < 600)[0]:
            out[i] = 2 / np.pi * np.sum(f * np.sin(t[i] * u))
        return out

    def N(self, t):
        s = np.zeros_like(t)
        for i in range(0, len(t), 4000):
            s[i:i + 4000] = -np.sin(np.outer(t[i:i + 4000], self.ln)) @ self.coef
        return theta(t) / np.pi + 1 + s / np.pi + self.pole(t) - 1


def crossings(t, N):
    Nm = np.maximum.accumulate(N)
    lo = max(1, int(np.ceil(N[0] + 0.5)))
    lev = np.arange(lo, int(np.floor(Nm[-1] - 0.5)) + 1) - 0.5
    i = np.clip(np.searchsorted(Nm, lev), 1, len(t) - 1)
    return t[i - 1] + (lev - N[i - 1]) / (N[i] - N[i - 1] + 1e-300) * (t[i] - t[i - 1])


def fit(ts, Ns, p0, U):
    x0 = np.append(p0, np.median(Ns - np.sum(kc(U * (ts[:, None] - p0[None, :])), axis=1)))

    def f(p):
        return p[-1] + np.sum(kc(U * (ts[:, None] - p[None, :-1])), axis=1) - Ns

    def J(p):
        return np.hstack([-U * k1(U * (ts[:, None] - p[None, :-1])), np.ones((len(ts), 1))])

    r = least_squares(f, x0, jac=J, method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=300)
    return np.sort(r.x[:-1])


def window(Z, T, X, nz=60, mg=25):
    P = Primes(X)
    k0 = int(np.searchsorted(Z, T))
    g = Z[k0:k0 + nz]
    d = 2 * np.pi / np.log(T / (2 * np.pi))
    t = np.arange(g[0] - (mg + 6) * d, g[-1] + (mg + 6) * d, min(d / 10, np.pi / (4 * P.U)))
    N = P.N(t)
    est = crossings(t, N)
    p0 = est[(est > g[0] - mg * d) & (est < g[-1] + mg * d)]
    sm = (t > g[0] - (mg - 5) * d) & (t < g[-1] + (mg - 5) * d)
    p = fit(t[sm], N[sm], p0, P.U)
    e0 = np.array([np.min(np.abs(est - z)) for z in g])
    e1 = np.array([np.min(np.abs(p - z)) for z in g])
    return e0, e1


def recover_all(X, TMAX):
    P = Primes(X)
    dmin = 2 * np.pi / np.log(TMAX / (2 * np.pi))
    t = np.arange(10.0, TMAX + 60, min(dmin / 10, np.pi / (4 * P.U)))
    N = P.N(t)
    est = crossings(t, N)
    out = np.full(len(est), np.nan)
    B, MG = 60, 25
    for b0 in range(0, len(est), B):
        a, b = max(0, b0 - MG), min(len(est), b0 + B + MG)
        p0 = est[a:b]
        d = 2 * np.pi / np.log(max(p0.mean(), 20) / (2 * np.pi))
        lo = 10.0 if a == 0 else p0[0] + 5 * d
        hi = p0[-1] + 3 * d if b == len(est) else p0[-1] - 5 * d
        sm = (t > lo) & (t < hi)
        p = fit(t[sm], N[sm], p0, P.U)
        m = min(B, len(est) - b0)
        out[b0:b0 + m] = p[b0 - a:b0 - a + m]
    return out, len(P.ln)


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'threshold':
        Z = np.fromfile(sys.argv[2], dtype=np.float64)
        for X in [30, 100, 300]:
            row = []
            for r in [2.0, 1.6, 1.3, 1.0]:
                T = 2 * np.pi * (X / r) ** 2
                if T > Z[-1] - 100:
                    continue
                e0, e1 = window(Z, T, X)
                row.append(f'X = {r} sqrt(T/2pi) (T = {T:.0f}): max error {e1.max():.1e} (crossings {e0.max():.1e})')
            print(f'X = {X}: ' + '; '.join(row), flush=True)
    else:
        X, TMAX, KEEP = float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
        Z = np.fromfile(sys.argv[5], dtype=np.float64)
        t0 = time.time()
        out, npp = recover_all(X, TMAX)
        keep = out[out <= KEEP]
        zt = Z[Z <= KEEP]
        print(f'prime powers up to {X:g}: {npp}; zeros recovered up to {KEEP:g}: {len(keep)} (true count {len(zt)}) '
              f'[{time.time() - t0:.0f} s]')
        if len(keep) == len(zt):
            e = np.abs(keep - zt)
            print(f'max error {e.max():.2e}, rms {np.sqrt(np.mean(e ** 2)):.2e}; above 1e-3: '
                  f'{[int(i) + 1 for i in np.nonzero(e > 1e-3)[0]]}')
        if len(sys.argv) > 6:
            keep.astype(np.float64).tofile(sys.argv[6])
