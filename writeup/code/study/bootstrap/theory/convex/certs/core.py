"""Core numerics for the certificate problem.

Closed form used everywhere (derived in README.md):
  Phi(u)   = 2 e^{u/2} sum_n (2y^2 - 3y) e^{-y},  y = pi n^2 e^{2u}
  (Phi*Phi)(u) = 4 sum_{N>=1} d(N) (pi N)^{-1/2} X^{1/2} R(X),   X = pi N e^u,
  R(X) = (4X^4 + 9X^2) K0(2X) - 12 X^3 K1(2X)
(each (m,n) pair of theta terms gives int z^{a-b-1} exp(-Az-C/z) dz = 2 (C/A)^{(a-b)/2} K_{a-b}(2 sqrt(AC));
 only N = mn enters, hence the divisor function d(N)).
d/du = X d/dX keeps the span {X^{1/2} P0(X) K0(2X) + X^{1/2} P1(X) K1(2X)}:
  D[X^{k+1/2} K0(2X)] = (k+1/2) X^{k+1/2} K0 - 2 X^{k+3/2} K1
  D[X^{k+1/2} K1(2X)] = (k-1/2) X^{k+1/2} K1 - 2 X^{k+3/2} K0
Basis: g_j = (-1)^j D^{2j} (Phi*Phi)  <->  ghat_j(t) = Xi(t)^2 t^{2j}.
"""
from fractions import Fraction
import mpmath as mp


def D(P0, P1):
    """apply X d/dX to X^{1/2}[P0 K0(2X) + P1 K1(2X)] (coefficient lists, Fractions)."""
    n = max(len(P0), len(P1)) + 1
    Q0 = [Fraction(0)] * n
    Q1 = [Fraction(0)] * n
    for k, a in enumerate(P0):
        if a:
            Q0[k] += a * (Fraction(k) + Fraction(1, 2))
            Q1[k + 1] += -2 * a
    for k, c in enumerate(P1):
        if c:
            Q1[k] += c * (Fraction(k) - Fraction(1, 2))
            Q0[k + 1] += -2 * c
    return Q0, Q1


def basis_polys(J, extra_derivs=0):
    """list over j=0..J of lists over r=0..extra_derivs of (P0,P1) for D^r g_j."""
    P0 = [Fraction(0), Fraction(0), Fraction(9), Fraction(0), Fraction(4)]
    P1 = [Fraction(0), Fraction(0), Fraction(0), Fraction(-12)]
    cur = (P0, P1)
    out = []
    for j in range(J + 1):
        sgn = -1 if j % 2 else 1
        lst = []
        tmp = cur
        for r in range(extra_derivs + 1):
            lst.append(([sgn * a for a in tmp[0]], [sgn * a for a in tmp[1]]))
            tmp = D(*tmp)
        out.append(lst)
        cur = D(*D(*cur))
    return out


def divisor_count(N):
    c = 0
    i = 1
    while i * i <= N:
        if N % i == 0:
            c += 1 if i * i == N else 2
        i += 1
    return c


class Basis:
    """Evaluate D^r g_j(u) for j=0..J, r=0..R at real u (mpmath, current mp.dps)."""

    def __init__(self, J, R=0):
        self.J, self.R = J, R
        self.fr = basis_polys(J, R)
        self.mp_polys = None
        self.dps = None

    def _prep(self):
        if self.dps != mp.mp.dps:
            self.mp_polys = [[([mp.mpf(a.numerator) / a.denominator for a in P0],
                               [mp.mpf(a.numerator) / a.denominator for a in P1])
                              for (P0, P1) in lst] for lst in self.fr]
            self.dps = mp.mp.dps
            self.maxdeg = max(max(len(P0), len(P1)) for lst in self.fr for (P0, P1) in lst)

    def eval(self, u, Nmax=None):
        """returns list[j][r] of D^r g_j(u)."""
        self._prep()
        u = mp.mpf(u)
        eu = mp.exp(u)
        out = [[mp.mpf(0)] * (self.R + 1) for _ in range(self.J + 1)]
        X1 = mp.pi * eu
        lim = (mp.mp.dps + 20) * mp.log(10)
        N = 1
        while True:
            X = mp.pi * N * eu
            if N >= 2 and Nmax is None:
                # crude bound of term N relative to term 1
                if 2 * (X - X1) - (self.maxdeg + 2) * mp.log(X / X1) - mp.log(N) > lim:
                    break
            if Nmax is not None and N > Nmax:
                break
            k0 = mp.besselk(0, 2 * X)
            k1 = mp.besselk(1, 2 * X)
            pref = 4 * divisor_count(N) / mp.sqrt(mp.pi * N) * mp.sqrt(X)
            pw = [mp.mpf(1)]
            for _ in range(self.maxdeg):
                pw.append(pw[-1] * X)
            for j in range(self.J + 1):
                for r in range(self.R + 1):
                    P0, P1 = self.mp_polys[j][r]
                    v = mp.fdot(P0, pw[:len(P0)]) * k0 + mp.fdot(P1, pw[:len(P1)]) * k1
                    out[j][r] += pref * v
            N += 1
        return out


# ---------------- Phi directly (for checks) -----------------
def Phi(u, nmax=None):
    u = mp.mpf(u)
    if u < 0:
        u = -u
    s = mp.mpf(0)
    n = 1
    while True:
        y = mp.pi * n * n * mp.exp(2 * u)
        t = (2 * y * y - 3 * y) * mp.exp(-y)
        s += t
        if y > (mp.mp.dps + 10) * 2.31 + 50 and n > 2:
            break
        n += 1
    return 2 * mp.exp(u / 2) * s


# ---------------- Xi side -----------------
def Xi(t):
    t = mp.mpf(t)
    return -(t * t + mp.mpf(1) / 4) / 2 * mp.pi ** (-mp.mpf(1) / 4) * abs(mp.gamma(mp.mpf(1) / 4 + 1j * t / 2)) * mp.siegelz(t)


def xi_moments(J, T=150, panel=1.0, nodes=24):
    """m_j = (1/2pi) int_R Xi^2 t^{2j} dt  and  a_j = (1/2pi) int_R Xi^2 t^{2j} (Re psi(1/4+it/2) - log pi) dt."""
    xs, ws = gauss_legendre(nodes)
    m = [mp.mpf(0)] * (J + 1)
    a = [mp.mpf(0)] * (J + 1)
    lp = mp.log(mp.pi)
    k = 0
    while k * panel < T:
        lo = k * panel
        for x, w in zip(xs, ws):
            t = lo + panel * (x + 1) / 2
            ww = w * panel / 2
            X2 = Xi(t) ** 2
            ps = mp.re(mp.digamma(mp.mpf(1) / 4 + 1j * t / 2)) - lp
            tp = mp.mpf(1)
            t2 = t * t
            for j in range(J + 1):
                m[j] += ww * X2 * tp
                a[j] += ww * X2 * tp * ps
                tp *= t2
        k += 1
    # (1/2pi) int_R = (1/pi) int_0^inf
    return [x / mp.pi for x in m], [x / mp.pi for x in a]


_gl_cache = {}


def gauss_legendre(n):
    key = (n, mp.mp.dps)
    if key not in _gl_cache:
        xs, ws = [], []
        # Newton on Legendre polynomials
        for i in range(1, n + 1):
            x = mp.cos(mp.pi * (i - mp.mpf(1) / 4) / (n + mp.mpf(1) / 2))
            for _ in range(100):
                p0, p1 = mp.mpf(1), x
                for k in range(2, n + 1):
                    p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
                dp = n * (x * p1 - p0) / (x * x - 1)
                dx = p1 / dp
                x -= dx
                if abs(dx) < mp.mpf(10) ** (-mp.mp.dps - 5):
                    break
            xs.append(x)
            ws.append(2 / ((1 - x * x) * dp * dp))
        _gl_cache[key] = (xs, ws)
    return _gl_cache[key]


def von_mangoldt(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n
            while m % p == 0:
                m //= p
            return mp.log(p) if m == 1 else mp.mpf(0)
    return mp.mpf(0)
