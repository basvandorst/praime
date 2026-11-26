"""Weil quadratic form D_P(f * f~) on L^2[-L/2, L/2], Legendre basis.

D_P = A - C_P with
  A(g) = 2 int g cosh(u/2) + a0 g(0) - int_0^inf k(u) [g(u)+g(-u)-2g(0)] du,
  k(u) = e^{-u/2}/(1-e^{-2u}),  a0 = psi(1/4) - log(pi)  (= Re psi(1/4) - log pi),
  C_P(g) = sum_{x in P, k>=1} log x x^{-k/2} (g(k log x) + g(-k log x)).
For g = f * f~, D(g) = int int f(x) conj f(y) D(x-y) dx dy.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss, legvander
from scipy.special import digamma

A0 = digamma(0.25) - np.log(np.pi)


def kfun(u):
    return np.exp(-u / 2) / (-np.expm1(-2 * u))


def ktail(L):
    # int_L^inf e^{-u/2}/(1-e^{-2u}) du = sum_j e^{-(2j+1/2)L}/(2j+1/2)
    j = np.arange(0, 2000)
    return np.sum(np.exp(-(2 * j + 0.5) * L) / (2 * j + 0.5))


def phi(n_max, L, x):
    """orthonormal Legendre basis on [-L/2, L/2]; returns array (len(x), n_max)."""
    V = legvander(2 * np.asarray(x) / L, n_max - 1)
    return V * np.sqrt((2 * np.arange(n_max) + 1) / L)


def corr(N, L, u, Q=None):
    """G[i, m, n] = g_mn(u_i) = int phi_m(x) phi_n(x - u_i) dx, 0 <= u_i <= L."""
    if Q is None:
        Q = N + 8
    t, w = leggauss(Q)
    u = np.atleast_1d(u)
    G = np.zeros((len(u), N, N))
    for i, ui in enumerate(u):
        lo, hi = ui - L / 2, L / 2
        if hi <= lo:
            continue
        x = (hi - lo) / 2 * t + (hi + lo) / 2
        wx = w * (hi - lo) / 2
        Pm = phi(N, L, x)
        Pn = phi(N, L, x - ui)
        G[i] = (Pm * wx[:, None]).T @ Pn
    return G


def form_A(N, L, Qu=None):
    """matrix of A on the basis."""
    if Qu is None:
        Qu = 2 * N + 60
    # pole part
    t, w = leggauss(N + 40)
    x = L / 2 * t
    wx = w * L / 2
    P = phi(N, L, x)
    al = (P * (wx * np.exp(x / 2))[:, None]).sum(0)
    be = (P * (wx * np.exp(-x / 2))[:, None]).sum(0)
    M = np.outer(al, be) + np.outer(be, al)
    # archimedean part
    tu, wu = leggauss(Qu)
    u = L / 2 * (tu + 1)
    wuu = wu * L / 2
    G = corr(N, L, u)
    S = G + np.transpose(G, (0, 2, 1)) - 2 * np.eye(N)[None]
    M += A0 * np.eye(N)
    M -= np.einsum('i,i,imn->mn', wuu, kfun(u), S)
    M += 2 * ktail(L) * np.eye(N)
    return M


def comb_atoms(P, L, kmax=200):
    """list of (position, weight) for comb atoms with 0 < position < L."""
    out = []
    for x in P:
        lx = np.log(x)
        for k in range(1, kmax + 1):
            if k * lx >= L:
                break
            out.append((k * lx, lx * x ** (-k / 2)))
    return out


def form_C(N, L, P):
    M = np.zeros((N, N))
    atoms = comb_atoms(P, L)
    if not atoms:
        return M
    pos = np.array([a[0] for a in atoms])
    wt = np.array([a[1] for a in atoms])
    G = corr(N, L, pos)
    M = np.einsum('i,imn->mn', wt, G + np.transpose(G, (0, 2, 1)))
    return M


def form_D(N, L, P):
    return form_A(N, L) - form_C(N, L, P)


def primes_upto(X):
    X = int(X)
    s = np.ones(X + 1, bool)
    s[:2] = False
    for p in range(2, int(X ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.nonzero(s)[0].astype(float)


def lam_min(M):
    return np.linalg.eigvalsh((M + M.T) / 2)[0]
