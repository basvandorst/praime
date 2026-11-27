"""Third reviewer's Weil-form engine (own code; imports nothing from the repo).

Explicit formula for Lambda(s) = q^{s/2} prod_j Gamma_R(s + kappa_j) L(s), -L'/L(s) = sum c(n) n^{-s} (c real here),
test function f on [-a, a] (u = log x), g(u) = int f(v+u) f(v) dv, h(t) = int g(u) e^{iut} du = |fhat(t)|^2:

  sum_gamma h(gamma) = [pole terms] + (log q) g(0) + sum_j A_{kappa_j}(g) - 2 sum_n c(n) n^{-1/2} g(log n)

  A_kappa(g) = -(gamma_E + log pi) g(0) + int_0^inf [2 e^{-2u} g(0) - 2 Re(e^{-(1/2+kappa) u}) g(u)] / (1 - e^{-2u}) du

(derivation: psi(w) = -gamma_E + int_0^inf (e^{-x} - e^{-wx})/(1 - e^{-x}) dx at w = (1/2 + kappa + it)/2, x = 2u,
psi_R(s) = -(1/2) log pi + (1/2) psi(s/2), integrated against h(t)/(2 pi); the dual side gives kappa-bar and g(-u) = g(u)
for real f.)  Since g = 0 for u >= 2a, the tail int_{2a}^inf 2 e^{-2u}/(1-e^{-2u}) du = -log(1 - e^{-4a}) is exact.
This route uses no digamma at all (the second reviewer's uses psi_R and a tail series): an independent evaluation.

Basis: orthonormal Legendre polynomials on [-a, a], degree < N; g(0) = |c|^2. Parity blocks: even / odd degree.
"""
import numpy as np
from numpy.polynomial import legendre as npl

EULER = 0.5772156649015328606065120900824


def gl(n, lo, hi):
    x, w = npl.leggauss(n)
    return 0.5 * (hi - lo) * x + 0.5 * (hi + lo), 0.5 * (hi - lo) * w


def basis(N, a, v):
    v = np.asarray(v, float)
    V = npl.legvander(np.clip(v / a, -1.0, 1.0), N - 1) * np.sqrt((2 * np.arange(N) + 1) / (2 * a))
    V[np.abs(v) > a * (1 + 1e-14)] = 0.0
    return V


def Smat(N, a, u):
    """S(u) = sym part of C(u), C_ij = int_{-a}^{a-u} phi_i(v+u) phi_j(v) dv (exact Gauss: polynomial integrand)."""
    if u >= 2 * a:
        return np.zeros((N, N))
    v, w = gl(N + 4, -a, a - u)
    A = basis(N, a, v + u)
    B = basis(N, a, v)
    C = (A * w[:, None]).T @ B
    return 0.5 * (C + C.T)


def primes(n):
    return [p for p in range(2, n + 1) if all(p % q for q in range(2, int(p ** 0.5) + 1))]


class Eng:
    def __init__(self, X, N, nu=700):
        self.X, self.N = float(X), int(N)
        self.a = 0.5 * np.log(X)
        self.u, self.w = gl(nu, 0.0, 2 * self.a)
        self.Su = np.array([Smat(self.N, self.a, x) for x in self.u])
        self.den = -np.expm1(-2 * self.u)
        self.I = np.eye(self.N)
        self._S = {}
        # prime powers p^m < X (p^m = X sits where g = 0)
        self.pp = [(p, m) for p in primes(int(X)) for m in range(1, 64) if p ** m < X * (1 - 1e-12)]
        ev = np.arange(0, self.N, 2)
        od = np.arange(1, self.N, 2)
        self.blocks = [ev, od]

    def S(self, u):
        k = round(float(u), 12)
        if k not in self._S:
            self._S[k] = Smat(self.N, self.a, u)
        return self._S[k]

    def arch(self, kappas):
        """sum_j A_{kappa_j} as an N x N real symmetric matrix (kappas: list of complex, closed under conj)."""
        M = np.zeros((self.N, self.N))
        for k in kappas:
            z = 0.5 + complex(k)
            wre = 2.0 * np.real(np.exp(-z * self.u))           # 2 Re e^{-(1/2+kappa)u}
            # integrand matrix: (2 e^{-2u} I - wre S(u)) / den ; integrate node-wise (combined, finite at u -> 0)
            Iu = (2 * np.exp(-2 * self.u) / self.den)[:, None, None] * self.I[None] - (wre / self.den)[:, None, None] * self.Su
            M += np.tensordot(self.w, Iu, axes=(0, 0))
            M += (-EULER - np.log(np.pi) - np.log1p(-np.exp(-4 * self.a))) * self.I
        return 0.5 * (M + M.T)

    def prime(self, p, m):
        """matrix of -2 log p p^{-m/2} g(m log p) (coefficient of c(p^m)/log p ... i.e. multiply by the weight w(p^m)/log p)."""
        return -2.0 * np.log(p) * p ** (-m / 2.0) * self.S(m * np.log(p))

    def center(self):
        v = np.zeros(self.N)
        v[0] = np.sqrt(2 * self.a)
        return v
