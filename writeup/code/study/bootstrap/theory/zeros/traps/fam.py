"""Fast high-precision evaluation of the family g = sum_j b_j g_j (python-flint arb, own K0/K1 from kbessel.py).

g_j = (-1)^j D^{2j}(Phi*Phi), D = d/du;  closed form (README of the certs folder):
  D^r g_j(u) = sum_N 4 d(N) (pi N)^{-1/2} X^{1/2} [P0_{jr}(X) K0(2X) + P1_{jr}(X) K1(2X)],  X = pi N e^u,
P0, P1 rational polynomials from core.basis_polys.

NUMERICS ONLY (not a proof): every input X is replaced by the exact midpoint of its ball before K0, K1 are
called (so the e^{2x} radius amplification of the K series branch does not eat the precision) and the N-sum
is truncated when the crude size e^{-2(X_N - X_1)} (X_N/X_1)^M N is below 2^-(prec+40).
"""
import math
from fractions import Fraction
import common  # noqa: F401  (sys.path)
from flint import arb, arb_poly, arb_mat, fmpq, ctx
from core import basis_polys, divisor_count
from kbessel import K01

ctx.prec = 900


def fq(x):
    x = Fraction(x)
    return arb(fmpq(x.numerator, x.denominator))


def mid(x):
    return arb(x.mid())


def astr(x, n=60):
    return x.mid().str(n, radius=False)


def a_from_str(s):
    return arb(s)


class Fam:
    def __init__(self, J, R):
        self.J, self.R = J, R
        self.fr = basis_polys(J, R)
        self.ap = [[(arb_poly([fq(a) for a in P0]), arb_poly([fq(a) for a in P1])) for (P0, P1) in lst]
                   for lst in self.fr]
        self.M = max(max(len(P0), len(P1)) for lst in self.fr for (P0, P1) in lst)

    def comb(self, b, r):
        """arb_poly pair of D^r (sum_j b_j g_j) (b: list of arb, len <= J+1)."""
        P0 = arb_poly([0])
        P1 = arb_poly([0])
        for j, bj in enumerate(b):
            P0 += self.ap[j][r][0] * bj
            P1 += self.ap[j][r][1] * bj
        return P0, P1

    def _terms(self, eu, Nmax=None):
        """yield (pref, X, k0, k1) for N = 1, 2, ... with eu = e^u (arb)."""
        PI = arb.pi()
        lim = (ctx.prec + 40) * math.log(2)
        X1 = None
        N = 1
        while True:
            X = mid(PI * N * eu)
            if N == 1:
                X1 = float(X)
            else:
                xf = float(X)
                if 2 * (xf - X1) - (self.M + 2) * math.log(xf / X1) - math.log(N) > lim:
                    return
            if Nmax is not None and N > Nmax:
                return
            k0, k1 = K01(2 * X)
            pref = 4 * divisor_count(N) / (PI * N).sqrt() * X.sqrt()
            yield pref, X, k0, k1
            N += 1

    def rows(self, eu, R=None, J=None):
        """list[j][r] of D^r g_j(u), eu = e^u as arb (exact n for u = log n)."""
        R = self.R if R is None else R
        J = self.J if J is None else J
        out = [[arb(0)] * (R + 1) for _ in range(J + 1)]
        for pref, X, k0, k1 in self._terms(eu):
            for j in range(J + 1):
                for r in range(R + 1):
                    P0, P1 = self.ap[j][r]
                    out[j][r] = out[j][r] + pref * (P0(X) * k0 + P1(X) * k1)
        return out

    def eval_polys(self, polys, eu):
        """values of the functions given by a list of (P0, P1) arb_poly pairs at e^u = eu."""
        out = [arb(0)] * len(polys)
        for pref, X, k0, k1 in self._terms(eu):
            for i, (P0, P1) in enumerate(polys):
                out[i] = out[i] + pref * (P0(X) * k0 + P1(X) * k1)
        return out


def dot(a, b):
    s = arb(0)
    for x, y in zip(a, b):
        s += x * y
    return s


def solve(rows, rhs):
    """solve sum_j A[i][j] x_j = rhs[i] (midpoints, numerics)."""
    n = len(rows)
    A = arb_mat([[mid(x) for x in r] for r in rows])
    B = arb_mat([[mid(x)] for x in rhs])
    X = A.solve(B)
    return [X[i, 0] for i in range(n)]
