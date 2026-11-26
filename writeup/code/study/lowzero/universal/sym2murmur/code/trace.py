"""Eichler-Selberg trace of T_n on S_k(Gamma0(N)), N an odd prime, gcd(n, N) = 1, k in {2, 4} (own code).
24*Tr T_n = 24*A1 - S2 + 24*A3 + 24*A4 with
  A1 = n^{k/2-1} (k-1)(N+1)/12 if n is a square,
  S2 = sum_{t^2 < 4n} U_{k-2}(t,n) * [12H(4n-t^2) (1 + chi_N(t^2-4n)) + N * 12H((4n-t^2)/N^2) [N^2 | 4n-t^2]],
  A3 = -sum_{d | n} min(d, n/d)^{k-1},   A4 = sigma_1(n) if k = 2.
U_0 = 1, U_2 = t^2 - n.  Then sum_f a_f(p)^2 = w^T G^{-1} w with G_{m m'} = sum_f a_f(m) a_f(m') = sum_{d|(m,m')} d^{k-1} t(m m'/d^2)
and w_m = sum_f a_f(p) a_f(m), over a set S of m coprime to N with G_S invertible (Hecke multiplicativity, newforms only)."""
from math import isqrt, gcd
import numpy as np
from sympy import divisors, divisor_sigma


def leg_table(N):
    leg = -np.ones(N, np.int64)
    leg[0] = 0
    leg[(np.arange(1, N, dtype=np.int64) ** 2) % N] = 1
    return leg


class Tracer:
    def __init__(self, N, k, H12):
        self.N, self.k, self.H = N, k, H12
        self.leg = leg_table(N)
        self.cache = {}

    def tr(self, n):
        if n in self.cache:
            return self.cache[n]
        N, k, H = self.N, self.k, self.H
        assert n % N != 0
        r = isqrt(4 * n)
        tmax = r if r * r < 4 * n else r - 1
        t = np.arange(-tmax, tmax + 1, dtype=np.int64)
        D = 4 * n - t * t
        term = H[D].astype(np.int64) * (1 + self.leg[(-D) % N])
        msk = (D % (N * N) == 0)
        if msk.any():
            term[msk] += N * H[D[msk] // (N * N)]
        if k == 2:
            S2 = int(term.sum())
        elif k == 4:
            S2 = int(((t * t - n) * term).sum())
        else:
            raise ValueError
        s = isqrt(n)
        A1x24 = 2 * s ** (k - 2) * (k - 1) * (N + 1) if s * s == n else 0
        A3 = -sum(min(d, n // d) ** (k - 1) for d in divisors(n))
        A4 = int(divisor_sigma(n)) if k == 2 else 0
        tot = A1x24 - S2 + 24 * (A3 + A4)
        assert tot % 24 == 0, (n, N, k, tot)
        self.cache[n] = tot // 24
        return tot // 24

    def pair(self, m1, m2):
        """sum_f a_f(m1) a_f(m2) (m1, m2 coprime to N)."""
        g = gcd(m1, m2)
        return sum(d ** (self.k - 1) * self.tr(m1 * m2 // (d * d)) for d in divisors(g))

    def basis(self, maxm=10 ** 5):
        dim = self.tr(1)
        S, G = [], np.zeros((0, 0))
        m = 0
        while len(S) < dim and m < maxm:
            m += 1
            if m % self.N == 0:
                continue
            S2 = S + [m]
            G2 = np.array([[self.pair(a, b) for b in S2] for a in S2], dtype=float)
            sv = np.linalg.svd(G2, compute_uv=False)
            if sv[-1] > 1e-9 * sv[0]:
                S, G = S2, G2
        assert len(S) == dim, (self.N, self.k, dim, S)
        self.S, self.G, self.dim = S, G, dim
        self.Gint = [[self.pair(a, b) for b in S] for a in S]
        return S

    def sum_ap2(self, p):
        """sum over newforms of a_f(p)^2 (p prime, p != N), exact (integrality checked)."""
        w = np.array([self.pair(p, m) for m in self.S], dtype=float)
        x = float(w @ np.linalg.solve(self.G, w))
        xi = round(x)
        return xi, abs(x - xi)
