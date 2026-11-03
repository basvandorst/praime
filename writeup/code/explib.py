import numpy as np
from math import log, sqrt, pi, exp

GAMMA = 0.5772156649015329
C2 = 0.6601618158468696
ZEROS_PATH = "/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/zeros_dl/zeros1"


def sieve(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return s


def primes_upto(n):
    return np.flatnonzero(sieve(n))


def segments(X, seg=10**8):
    base = primes_upto(int(sqrt(X + 8)) + 1)
    for lo in range(0, X, seg):
        hi = min(lo + seg, X)
        comp = np.zeros(hi - lo + 8, dtype=bool)
        for p in base:
            start = p * p - lo if lo == 0 else (-lo) % p
            if start < 0:
                start += ((-start) // p + 1) * p
            comp[start::p] = True
        if lo == 0:
            comp[:2] = True
        yield lo, ~comp


def zeros(n=100000):
    return np.loadtxt(ZEROS_PATH, max_rows=n)


def spectrum(values, weights, ts, chunk=20000):
    acc = np.zeros(len(ts), dtype=np.complex128)
    lv = np.log(values)
    w = weights / np.sqrt(values)
    for i in range(0, len(values), chunk):
        acc += (w[None, i:i + chunk] * np.exp(-1j * np.outer(ts, lv[i:i + chunk]))).sum(axis=1)
    return acc


def main_term(X, L, ts):
    return (X ** (0.5 - 1j * ts) - L ** (0.5 - 1j * ts)) / (0.5 - 1j * ts)


def psi_from_zeros(xs, gam, chunk=400):
    rho = 0.5 + 1j * gam
    out = np.empty(len(xs))
    for i in range(0, len(xs), chunk):
        x = xs[i:i + chunk]
        term = (np.sqrt(x)[:, None] * np.exp(1j * np.outer(np.log(x), gam)) / rho[None, :]).sum(axis=1)
        out[i:i + chunk] = x - 2 * term.real - log(2 * pi) - 0.5 * np.log1p(-1 / x**2)
    return out


def hl_constant(offsets, qmax=10**6):
    k = len(offsets)
    c = 1.0
    for q in primes_upto(qmax):
        nu = len(set(o % q for o in offsets))
        c *= (1 - nu / q) / (1 - 1 / q) ** k
    return c


def unfold(gam):
    return np.diff(gam) * np.log(gam[:-1] / (2 * pi)) / (2 * pi)


def smooth_N(t):
    return t / (2 * pi) * np.log(t / (2 * pi * np.e)) + 7 / 8


def lam_array(N):
    s = sieve(N)
    Lam = np.zeros(N + 1)
    for p in np.flatnonzero(s):
        q = p
        while q <= N:
            Lam[q] = log(p)
            q *= p
    return Lam
