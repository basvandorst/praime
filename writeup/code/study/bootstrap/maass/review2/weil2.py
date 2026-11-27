"""Second reviewer's own Weil-form engine (written from scratch; derivation in the docstring).

Explicit formula (derived by the contour shift, h(t) = int g(u) e^{iut} du, g(u) = int f(v+u) conj f(v) dv):
  sum_rho h(gamma_rho) = g(0) log q
        + sum_j [ 2 Re psi_R(1/2+k_j) g(0) + int_0^inf (e^{-(1/2+k_j)u}(g(0)-g(u)) + e^{-(1/2+conj k_j)u}(g(0)-g(-u)))/(1-e^{-2u}) du ]
        - sum_n n^{-1/2} [ c(n) g(log n) + conj c(n) g(-log n) ],          -L'/L = sum c(n) n^{-s},
  Lambda(s) = q^{s/2} prod_j Gamma_R(s+k_j) L(s), Gamma_R(s) = pi^{-s/2} Gamma(s/2), psi_R = Gamma_R'/Gamma_R.
(used: psi(z) = -gamma + int_0^inf (e^{-x} - e^{-zx})/(1-e^{-x}) dx with z = (1/2+k+it)/2, x = 2u.)

Basis: orthonormal Legendre phi_i(v) = sqrt((2i+1)/(2a)) P_i(v/a) on [-a, a], a = log(X)/2.
Matrix convention: g(u) = c^H B(u) c with B(u) = C(u)^T, C(u)_ij = int phi_i(v+u) phi_j(v) dv.
Real symmetric problems use S(u) = (C(u)+C(u)^T)/2 and split into even/odd Legendre blocks.
"""
import numpy as np
from numpy.polynomial import legendre as L
import mpmath as mp

EG = 0.57721566490153286060651209


def gauss(n, lo, hi):
    x, w = L.leggauss(n)
    return 0.5 * (hi - lo) * x + 0.5 * (hi + lo), 0.5 * (hi - lo) * w


def phis(N, a, v):
    v = np.atleast_1d(np.asarray(v, dtype=float))
    out = np.zeros((len(v), N))
    inside = np.abs(v) <= a + 1e-15
    V = L.legvander(np.clip(v / a, -1, 1), N - 1) * np.sqrt((2 * np.arange(N) + 1) / (2 * a))
    out[inside] = V[inside]
    return out


def Cmat(N, a, u):
    """C(u)_ij = int phi_i(v+u) phi_j(v) dv for 0 <= u <= 2a (polynomial integrand: Gauss exact)."""
    if u >= 2 * a:
        return np.zeros((N, N))
    v, w = gauss(N + 3, -a, a - u)
    A = phis(N, a, v + u)
    B = phis(N, a, v)
    return (A * w[:, None]).T @ B


def psiR(z):
    return -0.5 * mp.log(mp.pi) + 0.5 * mp.digamma(z / 2)


def tail(z, a):
    """int_{2a}^inf e^{-z u}/(1-e^{-2u}) du = sum_m e^{-(z+2m)2a}/(z+2m)."""
    s = mp.mpc(0)
    m = 0
    while True:
        term = mp.exp(-(z + 2 * m) * 2 * a) / (z + 2 * m)
        s += term
        if abs(term) < 1e-25 or m > 2000:
            return s
        m += 1


class Engine:
    def __init__(self, X, N, nu=800):
        self.X, self.N = float(X), N
        self.a = 0.5 * np.log(X)
        self.u, self.wu = gauss(nu, 0.0, 2 * self.a)
        self.Cu = np.array([Cmat(N, self.a, uu) for uu in self.u])   # (nu, N, N)
        self.I = np.eye(N)
        self.cache = {}

    def C(self, u):
        k = round(float(u), 13)
        if k not in self.cache:
            self.cache[k] = Cmat(self.N, self.a, u)
        return self.cache[k]

    def S(self, u):
        c = self.C(u)
        return 0.5 * (c + c.T)

    def arch_herm(self, kappas):
        """Hermitian matrix (c^H M c) of sum_j A_{k_j}; g(u) <-> C(u)^T, g(-u) <-> C(u)."""
        N = self.N
        M = np.zeros((N, N), dtype=complex)
        den = -np.expm1(-2 * self.u)
        for k in kappas:
            z = mp.mpf(0.5) + mp.mpc(k)
            const = 2 * mp.re(psiR(z)) + 2 * mp.re(tail(z, self.a))   # u > 2a: g = 0, so g(0)-g = g(0)
            e1 = np.exp(-complex(z) * self.u) / den
            e2 = np.exp(-complex(z).conjugate() * self.u) / den
            IC = self.I[None] - self.Cu
            ICT = self.I[None] - np.transpose(self.Cu, (0, 2, 1))
            M += float(const) * self.I
            M += np.tensordot(self.wu * e1, ICT, axes=(0, 0)) + np.tensordot(self.wu * e2, IC, axes=(0, 0))
        return M

    def arch_sym(self, kappas):
        M = self.arch_herm(kappas)
        Ms = 0.5 * (M + M.conj().T)
        assert np.max(np.abs(Ms.imag)) < 1e-9 * max(1, np.max(np.abs(Ms))), "not real"
        return Ms.real

    def prime_sym(self, n, c):
        """real c(n): matrix of -n^{-1/2} [c g(log n) + c g(-log n)] = -2 c n^{-1/2} S(log n)."""
        return -2.0 * c * n ** -0.5 * self.S(np.log(n))

    def prime_herm(self, n, c):
        Cn = self.C(np.log(n))
        return -(n ** -0.5) * (c * Cn.T + np.conj(c) * Cn)

    def blocks(self):
        return [np.arange(0, self.N, 2), np.arange(1, self.N, 2)]

    def center_vec(self):
        """int phi_i over [-a,a]: only P_0 contributes, sqrt(2a)."""
        v = np.zeros(self.N)
        v[0] = np.sqrt(2 * self.a)
        return v


def primes_upto(X):
    ps = []
    for p in range(2, int(X) + 1):
        if all(p % q for q in range(2, int(p ** 0.5) + 1)):
            ps.append(p)
    return ps


def prime_powers(X):
    """(p, m) with p^m < X (p^m = X sits where g vanishes)."""
    out = []
    for p in primes_upto(X):
        m = 1
        while p ** m < X - 1e-12:
            out.append((p, m))
            m += 1
    return out
