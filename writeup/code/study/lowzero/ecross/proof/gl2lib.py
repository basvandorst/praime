"""gl2lib: zeros <-> primes for the one-point sum of an elliptic curve L-function (own code, written for the proof note).

Normalisation: L(s,E) = sum lambda(n) n^-s, lambda(n) = a_n/sqrt(n); Lambda(s) = N^{s/2} Gamma_C(s+1/2) L(s,E) = w Lambda(1-s);
X(s) = N^{1/2-s} Gamma_C(3/2-s)/Gamma_C(s+1/2), L(s) = w X(s) L(1-s);  X(1/2+it)^{-1} = e^{2 i theta(t)},
theta(t) = (t/2) ln N - t ln 2pi + Im lnGamma(1+it).
Explicit formula (proof note, Prop. 1): sum_rho h(rho) X(rho)^{-1} = -(1/2pi) sum_n Lambda_E(n) n^{-1/2} I_n + E_0,
I_n = int H(t) e^{i(2 theta(t) - t ln n)} dt.  Stationary phase at 2 theta'(t_n) = ln n.
"""
import math, json, os
import numpy as np
from scipy.special import loggamma, digamma, erf

TWO_PI = 2 * math.pi
REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', '..', '..'))
ZEROS = os.path.join(REPO, 'writeup/data/study/r144/fresh_zeros.npz')
CURVES = os.path.join(REPO, 'writeup/code/study/lowzero/deg2/fresh_set.json')


# ----------------------------------------------------------------------------- theta and its derivatives
def theta(N, t):
    t = np.asarray(t, float)
    return t / 2 * math.log(N) - t * math.log(TWO_PI) + np.imag(loggamma(1 + 1j * t))


def Lfreq(N, t):
    """2 theta'(t) = ln N - 2 ln 2pi + 2 Re psi(1+it)."""
    t = np.asarray(t, float)
    return math.log(N) - 2 * math.log(TWO_PI) + 2 * np.real(digamma(1 + 1j * t))


def polygamma_c(k, z, shift=12):
    """psi^{(k)}(z), k = 1, 2, 3, complex z with Re z >= 1: recurrence z -> z + shift, then the asymptotic series."""
    z = np.asarray(z, complex)
    acc = 0
    fk = math.factorial(k)
    for j in range(shift):
        acc = acc + (-1) ** (k + 1) * fk / (z + j) ** (k + 1)
    y = z + shift
    # psi^{(k)}(y) ~ (-1)^{k+1} [ (k-1)!/y^k + k!/(2 y^{k+1}) + sum_j B_{2j} (2j+k-1)!/((2j)! y^{2j+k}) ]
    B = [1 / 6, -1 / 30, 1 / 42, -1 / 30, 5 / 66]
    s = math.factorial(k - 1) / y ** k + fk / (2 * y ** (k + 1))
    for j, b in enumerate(B, start=1):
        s = s + b * math.factorial(2 * j + k - 1) / math.factorial(2 * j) / y ** (2 * j + k)
    return acc + (-1) ** (k + 1) * s


def theta_derivs(t):
    """2 theta'', 2 theta''', 2 theta'''' at real t (independent of N)."""
    z = 1 + 1j * np.asarray(t, float)
    f2 = -2 * np.imag(polygamma_c(1, z))          # d/dt 2 Re psi(1+it) = -2 Im psi'(1+it)
    f3 = -2 * np.real(polygamma_c(2, z))          # d/dt of that
    f4 = 2 * np.imag(polygamma_c(3, z))
    return f2, f3, f4


def stationary_points(N, lnm, t_guess=None, iters=6):
    """t with 2 theta'(t) = lnm (Newton)."""
    t = np.sqrt(np.exp(lnm) / N) * TWO_PI if t_guess is None else t_guess
    for _ in range(iters):
        f = Lfreq(N, t) - lnm
        f2, _, _ = theta_derivs(t)
        t = t - f / f2
    return t


# ----------------------------------------------------------------------------- weights
class Weight:
    """H(t) = (W * g_Delta)(t), W = indicator of [T1,T2] ('box') or Hann on [T1,T2] ('hann'); Delta = 0 means no smoothing.
    Provides H, H', H'' on the real line."""

    def __init__(self, kind, T1, T2, Delta=0.0):
        self.kind, self.T1, self.T2, self.D = kind, float(T1), float(T2), float(Delta)

    def support(self):
        m = 6 * self.D
        return self.T1 - m, self.T2 + m

    def __call__(self, t, d=0):
        t = np.asarray(t, float)
        T1, T2, D = self.T1, self.T2, self.D
        if self.kind == 'box':
            if D == 0:
                if d == 0:
                    return np.where((t > T1) & (t <= T2), 1.0, 0.0)
                return np.zeros_like(t)
            g = lambda x: np.exp(-(x / D) ** 2) / (math.sqrt(math.pi) * D)
            if d == 0:
                return 0.5 * (erf((t - T1) / D) - erf((t - T2) / D))
            if d == 1:
                return g(t - T1) - g(t - T2)
            if d == 2:
                gp = lambda x: -2 * x / D ** 2 * g(x)
                return gp(t - T1) - gp(t - T2)
        if self.kind == 'hann':
            L = T2 - T1; a = math.pi / L
            inside = (t > T1) & (t < T2)
            if D == 0:
                if d == 0: return np.where(inside, np.sin(a * (t - T1)) ** 2, 0.0)
                if d == 1: return np.where(inside, a * np.sin(2 * a * (t - T1)), 0.0)
                if d == 2: return np.where(inside, 2 * a * a * np.cos(2 * a * (t - T1)), 0.0)
            # Gaussian-smoothed Hann: W = (1 - cos(2a(t-T1)))/2 on [T1,T2]; convolution in closed form via erf with complex args
            # is avoided: numerical convolution on a fine grid (used only in the checks, small D)
            raise NotImplementedError
        raise ValueError(self.kind)

    def integral(self):
        if self.kind == 'box':
            return self.T2 - self.T1
        return (self.T2 - self.T1) / 2


# ----------------------------------------------------------------------------- curves and coefficients
def load_curves():
    d = json.load(open(CURVES))
    return {c['id']: c for c in d['curves']}


def load_zeros(cid):
    z = np.load(ZEROS)
    return np.sort(np.asarray(z[cid], float))


def prime_power_data(pari, ainv, bad, X):
    """n = p^k <= X (sorted) and Lambda_E(n) = (alpha^k + beta^k) ln p (good p), lambda_p^k ln p (bad p, multiplicative)."""
    pari('E__ = ellinit(%s); P__ = primes([2, %d]);' % (str(list(ainv)), int(X)))
    ap = np.array(pari('Vecsmall(apply(p -> ellap(E__, p), P__))'), dtype=float)
    pr = np.array(pari('Vecsmall(P__)'), dtype=float)
    lam = ap / np.sqrt(pr)
    isbad = np.isin(pr.astype(np.int64), np.array(bad, np.int64))
    # Chebyshev recursion for good p: c_k = alpha^k + beta^k, c_0 = 2, c_1 = lam, c_{k+1} = lam c_k - c_{k-1}
    ns, Ls = [], []
    c_prev, c_cur = np.full_like(lam, 2.0), lam.copy()
    k = 1
    while True:
        m = pr ** k <= X
        if not m.any():
            break
        Lam = np.where(isbad, lam ** k, c_cur) * np.log(pr)
        ns.append(pr[m] ** k); Ls.append(Lam[m])
        c_prev, c_cur = c_cur, lam * c_cur - c_prev
        k += 1
    n = np.concatenate(ns); L = np.concatenate(Ls)
    o = np.argsort(n)
    return n[o], L[o], dict(zip(pr.astype(np.int64)[:50].tolist(), ap[:50].astype(int).tolist()))


def lambda_all(pari, ainv, X):
    """lambda(n) = a_n/sqrt(n) for n = 1..X (PARI ellan)."""
    pari('E__ = ellinit(%s);' % str(list(ainv)))
    a = np.array(pari('Vecsmall(ellan(E__, %d))' % int(X)), dtype=float)
    n = np.arange(1, len(a) + 1, dtype=float)
    return n, a / np.sqrt(n)


# ----------------------------------------------------------------------------- the two sides
def zero_side(N, gam, H, r=1.0):
    gam = np.asarray(gam, float)
    return np.sum(H(gam) * np.exp(2j * theta(N, gam)) * np.exp(-1j * gam * math.log(r)))


def prime_side(N, n, coef, H, r=1.0, mode='exact', second_order=False, chunk=150000):
    """chunked wrapper of prime_side_1"""
    n = np.asarray(n, float); coef = np.asarray(coef, float)
    tot = 0j
    for s in range(0, len(n), chunk):
        tot += prime_side_1(N, n[s:s + chunk], coef[s:s + chunk], H, r, mode, second_order)
    return tot


def prime_side_1(N, n, coef, H, r=1.0, mode='exact', second_order=False):
    """-(1/2pi) sum coef(n) n^{-1/2} [stationary-phase value of int H e^{i(2theta - t ln(nr))}].
    mode 'exact': exact theta at the exact stationary point; 'simple': -(e^{3pi i/4}/sqrt2) N^{-1/4} r^{1/4}
    coef n^{-1/4} e(-2 sqrt(nr/N)) H(2 pi sqrt(nr/N))."""
    n = np.asarray(n, float)
    lo, hi = H.support()
    tg = TWO_PI * np.sqrt(n * r / N)
    m = (tg > lo - 2) & (tg < hi + 2)
    n, coef, tg = n[m], coef[m], tg[m]
    if mode == 'simple':
        tt = tg
        val = -(np.exp(3j * math.pi / 4) / math.sqrt(2)) * N ** -0.25 * r ** 0.25 * coef * n ** -0.25 \
            * np.exp(-2j * math.pi * 2 * np.sqrt(n * r / N)) * H(tt)
        return np.sum(val)
    lnm = np.log(n * r)
    t0 = stationary_points(N, lnm, tg)
    f2, f3, f4 = theta_derivs(t0)
    ph = 2 * theta(N, t0) - t0 * lnm
    g = H(t0)
    amp = g
    if second_order:
        g1, g2 = H(t0, 1), H(t0, 2)
        amp = g + (1j / (2 * f2)) * (g2 - g1 * f3 / f2 - g * f4 / (4 * f2) + 5 * g * f3 ** 2 / (12 * f2 ** 2))
    val = -(1 / TWO_PI) * coef * n ** -0.5 * np.sqrt(TWO_PI / f2) * np.exp(1j * (ph + math.pi / 4)) * amp
    return np.sum(val)
