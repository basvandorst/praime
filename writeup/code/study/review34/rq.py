"""Review 34: own code for the round-100 swing statistic. Nothing imported from the study scripts."""
import numpy as np, math
from scipy.linalg import eigh_tridiagonal

ZF = '/home/bas/praime/writeup/data/riemann/zeros_10000_flint.txt'
Z = np.array([float(l) for l in open(ZF) if l.strip()])


def Nbar(t):
    t = np.asarray(t, float)
    return t / (2 * math.pi) * np.log(t / (2 * math.pi * math.e)) + 7 / 8


def Ninv(y):
    """inverse of Nbar on t >= 2 pi by Newton"""
    y = np.asarray(y, float)
    t = np.maximum(2 * math.pi * math.e * np.ones_like(y), 2 * math.pi * (y + 1) / np.maximum(np.log(np.maximum(y, 2) + 1), 1))
    t = np.maximum(t, 2 * math.pi * 1.0001)
    for _ in range(60):
        f = Nbar(t) - y; fp = np.log(t / (2 * math.pi)) / (2 * math.pi)
        t = np.maximum(t - f / np.maximum(fp, 1e-6), 2 * math.pi * 1.000001)
    return t


def c2(g, d):
    rho = 0.5 + 1j * g
    return np.abs(((1 + d) ** rho - 1) / rho) ** 2


def rows(g, ds, chunk=250):
    """for positive points g (sorted), return dict d -> w_j with Q = 2 sum_j w_j,
    w_j = a_j sum_{k != j} a_k L(g_j - g_k) + a_j sum_k a_k L(g_j + g_k), L = 1/(4 + x^2)"""
    g = np.asarray(g, float); n = len(g)
    A = {d: c2(g, d) for d in ds}
    W = {d: np.zeros(n) for d in ds}
    for i0 in range(0, n, chunk):
        i1 = min(n, i0 + chunk)
        dm = g[i0:i1, None] - g[None, :]
        Lm = 1.0 / (4 + dm * dm)
        ii = np.arange(i0, i1); Lm[ii - i0, ii] = 0.0
        dp = g[i0:i1, None] + g[None, :]
        Lp = 1.0 / (4 + dp * dp)
        for d in ds:
            a = A[d]
            W[d][i0:i1] = a[i0:i1] * (Lm @ a + Lp @ a)
        del dm, Lm, dp, Lp
    return W


def rms(g, ds):
    W = rows(g, ds)
    return {d: 2 / d * math.sqrt(2 * 2 * W[d].sum()) for d in ds}


def gue_unf(seed, n=40000):
    """beta = 2 Hermite (Dumitriu-Edelman), H = tridiag(N(0,2), chi_{2k}) / sqrt 2: semicircle radius 2 sqrt n (exact)."""
    rng = np.random.default_rng(seed)
    dg = rng.normal(0, math.sqrt(2), n) / math.sqrt(2)
    off = np.sqrt(rng.chisquare(2 * np.arange(n - 1, 0, -1))) / math.sqrt(2)
    ev = eigh_tridiagonal(dg, off, eigvals_only=True)
    x = np.clip(np.sort(ev) / (2 * math.sqrt(n)), -1, 1)
    return n * (0.5 + (x * np.sqrt(1 - x * x) + np.arcsin(x)) / math.pi)


def gue_set(seed, lo, hi, n=40000):
    """stationary: random origin in the central bulk, points with lo + v in [lo, hi] in mean-count space, mapped by Ninv"""
    u = gue_unf(seed, n)
    rng = np.random.default_rng(10_000 + seed)
    o = rng.uniform(0.30 * n, 0.35 * n)
    v = u - o; v = v[(v >= 0) & (v <= hi - lo)]
    return Ninv(lo + v)
