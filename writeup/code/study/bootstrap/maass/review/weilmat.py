"""Own implementation (reviewer) of truncated Weil quadratic forms, x-space (Weil's form).

Test function phi in L^2[-a,a], orthonormal Legendre basis phi_i(x) = sqrt((2i+1)/(2a)) P_i(x/a).
g = phi * phi~ (autocorrelation), supported in [-2a, 2a], h = |phi^|^2.
For Lambda(s) = N^{s/2} prod_j Gamma_R(s+mu_j) L(s), -L'/L = sum Lambda_f(n) n^-s, Re(1/2+mu_j) > 0:

  sum_rho h(gamma) = log N g(0) + A[g] - 2 sum_n Lambda_f(n) n^{-1/2} g(log n)      (no pole)
  A[g] = g(0) sum_j 2 Re psi_R(1/2+mu_j) + int_0^inf K(u) (g(0) - g(u)) du
  K(u) = sum_j 2 Re exp(-(1/2+mu_j) u) / (1 - exp(-2u)),  psi_R(s) = -log(pi)/2 + psi(s/2)/2

(derived from psi(z) = int_0^inf [e^-x/x - e^{-zx}/(1-e^-x)] dx, x = 2u).
Matrix of the form: M = c0 I + int_0^{2a} K(u) (I - G(u)) du - sum_n 2 Lambda_f(n) n^{-1/2} G(log n),
c0 = sum_j 2 Re psi_R(1/2+mu_j) + int_{2a}^inf K, G(u)_ij = symmetrised int phi_i(y) phi_j(y-u) dy.
"""
import numpy as np
from numpy.polynomial import legendre as npl
import mpmath as mp


def gl(n, lo, hi):
    x, w = npl.leggauss(n)
    return 0.5 * (hi - lo) * x + 0.5 * (hi + lo), 0.5 * (hi - lo) * w


def basis_vals(N, a, x):
    """phi_i(x) for i < N, shape (len(x), N)."""
    V = npl.legvander(np.asarray(x) / a, N - 1)
    return V * np.sqrt((2 * np.arange(N) + 1) / (2 * a))


def Gmat(N, a, u, ny=None):
    """Symmetrised correlation matrix G(u), u in [0, 2a]."""
    if ny is None:
        ny = N + 4
    y, w = gl(ny, u - a, a)
    P1 = basis_vals(N, a, y)
    P2 = basis_vals(N, a, y - u)
    g = (P1 * w[:, None]).T @ P2
    return 0.5 * (g + g.T)


class Geo:
    """Geometry for a support: a = log(X)/2, basis size N, u-quadrature."""

    def __init__(self, X, N, nu=None):
        self.X, self.N = X, N
        self.a = 0.5 * np.log(X)
        if nu is None:
            nu = max(500, 4 * N)
        self.u, self.wu = gl(nu, 0.0, 2 * self.a)
        self.IG = np.array([np.eye(N) - Gmat(N, self.a, uu) for uu in self.u])  # (nu, N, N)
        self.Gcache = {}

    def G(self, u):
        key = round(u, 14)
        if key not in self.Gcache:
            self.Gcache[key] = Gmat(self.N, self.a, u)
        return self.Gcache[key]


def psiR(s):
    return -0.5 * mp.log(mp.pi) + 0.5 * mp.digamma(s / 2)


def c0_of(mus, a):
    s = mp.mpf(0)
    for mu in mus:
        mu = mp.mpc(mu)
        s += 2 * mp.re(psiR(mp.mpf(0.5) + mu))
        # tail int_{2a}^inf 2 Re e^{-(1/2+mu)u}/(1-e^{-2u}) du
        m = 0
        while True:
            z = mp.mpf(0.5) + mu + 2 * m
            term = 2 * mp.re(mp.exp(-z * 2 * a) / z)
            s += term
            if abs(term) < mp.mpf(10) ** (-30) or m > 400:
                break
            m += 1
    return float(s)


def Kfun(mus, u):
    num = np.zeros_like(u)
    for mu in mus:
        mu = complex(mu)
        num += 2 * np.exp(-(0.5 + mu.real) * u) * np.cos(mu.imag * u)
    return num / (-np.expm1(-2 * u))


def arch(geo, mus):
    K = Kfun(mus, geo.u)
    M = np.tensordot(geo.wu * K, geo.IG, axes=(0, 0))
    return c0_of(mus, geo.a) * np.eye(geo.N) + M


def prime_powers(X):
    out = []
    for p in range(2, int(X) + 1):
        if all(p % q for q in range(2, int(p ** 0.5) + 1)):
            m = 1
            while p ** m < X:   # p^m = X sits at the boundary, g(2a) = 0
                out.append((p, m))
                m += 1
    return out


def gammas(kind, r, eps):
    if kind == 1:
        return [eps + 1j * r, eps - 1j * r]
    if kind == 2:
        return [0, 2j * r, -2j * r]
    if kind == 3:
        return [eps + 3j * r, eps - 3j * r, eps + 1j * r, eps - 1j * r]
    if kind == 4:
        return [0, 2j * r, -2j * r, 4j * r, -4j * r]
    raise ValueError


def sym_chars(k, m):
    """Sym^k coefficient at p^m: sum_{j=0}^k cos((k-2j) m theta) -> dict l -> multiplicity."""
    d = {}
    for j in range(k + 1):
        l = abs(k - 2 * j) * m
        d[l] = d.get(l, 0) + 1
    return d


def form_matrix(geo, k, r, eps, moments, logN=0.0):
    """moments: dict p -> array c with c[l] = E cos(l theta), c[0] = 1."""
    M = arch(geo, gammas(k, r, eps)) + logN * np.eye(geo.N)
    for p, m in prime_powers(geo.X):
        w = 2 * np.log(p) * p ** (-m / 2)
        coef = sum(mult * moments[p][l] for l, mult in sym_chars(k, m).items())
        M -= w * coef * geo.G(m * np.log(p))
    return M


def true_moments(ap, L=40):
    mom = {}
    for p, a in ap.items():
        th = np.arccos(np.clip(a / 2, -1, 1))
        mom[int(p)] = np.cos(np.arange(L + 1) * th)
    return mom
