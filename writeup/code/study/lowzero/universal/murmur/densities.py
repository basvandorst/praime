"""Murmuration densities (weight k = 2) and the family-theorem kernel. Imported by the check scripts
(own code; nothing from the repository).

M_harm(y)  : harmonic (Petersson) density, prime level:      2 pi sqrt(y) sum_m mu^2(m)/(m phi(m)) J_1(4 pi sqrt(y)/m)
M_zub(y)   : Zubrilina Thm 1.1, k = 2, squarefree levels:     alpha sum_{1<=r<=2 sqrt y} nu(r) sqrt(4y-r^2) + beta sqrt y - gamma y
M_prime(y) : natural weights, prime level, DERIVED here from the one-point law with lines + Plancherel:
             (12/pi) sqrt(y) sum_{d sqfree} (w(d)/d) sum_{e|d} (mu(e)/e) F(d/e),  w(d) = prod_{p|d} p^2/((p-1)^2 (p+1)),
             F(delta) = sum_{|n| <= 2 sqrt(y)/delta} sqrt(1 - n^2 delta^2/(4y)) - pi sqrt(y)/delta   (Poisson form of
             sum_s J_1(4 pi s sqrt(y)/delta)/s);  for d > 2 sqrt y the inner sum is phi(d)/d.
G(y)       : kernel of the family theorem after the change of variables u = N y (N cancels):
             G(y) = int H(t) exp(i[2 Im lnGamma(1+it) - 2 t ln(2 pi) - t ln y]) dt,
             family value  E[eps Z(H)] ~ -(1/2pi) int_0^inf y^{-1} M(y) G(y) dy.
"""
import numpy as np
from scipy.special import j1, loggamma
import sympy as sp

PR = list(sp.primerange(2, 3_000_000))


def euler(f):
    s = 0.0
    for q in PR:
        s += np.log(f(float(q)))
    return float(np.exp(s))


ALPHA = 2 * np.pi * euler(lambda q: (1 - q - 2 * q ** 2 + q ** 4) / (q ** 4 - 2 * q ** 2 + q))
BETA = 2 * np.pi * euler(lambda q: (-1 + q ** 2 + q ** 3) / (q * (-1 + q + q ** 2)))
GAMMA = 12 * euler(lambda q: q * (1 + q) / (-1 + q + q ** 2))


def factor_sqfree(n):
    f = sp.factorint(n)
    if any(e > 1 for e in f.values()):
        return None
    return list(f)


# ---------- Zubrilina, k = 2 ----------
_NU = {}


def nu(r):
    if r not in _NU:
        v = 1.0
        for q in sp.factorint(r):
            v *= 1 + q ** 2 / (q ** 4 - 2 * q ** 2 - q + 1)
        _NU[r] = v
    return _NU[r]


def M_zub(y):
    y = np.asarray(y, dtype=float)
    out = BETA * np.sqrt(y) - GAMMA * y
    rmax = int(np.floor(2 * np.sqrt(y.max()))) + 1
    for r in range(1, rmax + 1):
        arg = 4 * y - r * r
        out = out + ALPHA * nu(r) * np.sqrt(np.clip(arg, 0, None))
    return out


# ---------- harmonic, prime level ----------
def M_harm(y, mmax=20000):
    y = np.asarray(y, dtype=float)
    out = np.zeros_like(y)
    sy = np.sqrt(y)
    for m in range(1, mmax + 1):
        f = factor_sqfree(m)
        if f is None:
            continue
        phi = 1
        for q in f:
            phi *= q - 1
        out += j1(4 * np.pi * sy / m) / (m * phi)
    # tail m > mmax: J_1(x) ~ x/2 ; sum_{m>mmax} mu^2/(m phi(m)) * 2 pi sqrt(y)/m  ~ tiny, ignored
    return 2 * np.pi * sy * out


# ---------- natural weights, prime level (derived) ----------
def _w(f):
    v = 1.0
    for q in f:
        v *= q * q / ((q - 1) ** 2 * (q + 1))
    return v


def _F(delta, y):
    sy = np.sqrt(y)
    nmax = int(np.floor(2 * sy.max() / delta)) + 1
    out = -np.pi * sy / delta
    for n in range(-nmax, nmax + 1):
        a = 1 - (n * delta) ** 2 / (4 * y)
        out = out + np.sqrt(np.clip(a, 0, None))
    return out


ZETA2 = np.pi ** 2 / 6
# c(delta) = (w(delta)/delta) * C0 / prod_{p|delta} (1 - u_p),  u_p = 1/((p-1)^2 (p+1)),  C0 = prod_p (1 - u_p)
C0 = euler(lambda q: 1 - 1 / ((q - 1) ** 2 * (q + 1)))


def c_delta(f, delta):
    v = C0 / delta
    for q in f:
        u = 1 / ((q - 1) ** 2 * (q + 1))
        v *= (q * q / ((q - 1) ** 2 * (q + 1))) / (1 - u)
    return v


_NUP = {}


def nu_prime(r):
    """nu'(r) = sum_{delta | r, delta squarefree} c(delta)."""
    if r not in _NUP:
        fr = list(sp.factorint(r))
        tot = 0.0
        for k in range(1 << len(fr)):
            sub = [q for i, q in enumerate(fr) if k >> i & 1]
            delta = 1
            for q in sub:
                delta *= q
            tot += c_delta(sub, delta)
        _NUP[r] = tot
    return _NUP[r]


def M_prime(y):
    """Closed form (Poisson in s, regrouped by delta = d/e and r = n delta):
       M_prime(y) = 2 pi sqrt(y) - 12 y + (12/pi) sum_{1 <= r <= 2 sqrt y} nu'(r) sqrt(4y - r^2)."""
    y = np.asarray(y, dtype=float)
    out = 2 * np.pi * np.sqrt(y) - 12 * y
    rmax = int(np.floor(2 * np.sqrt(y.max()))) + 1
    for r in range(1, rmax + 1):
        out = out + (12 / np.pi) * nu_prime(r) * np.sqrt(np.clip(4 * y - r * r, 0, None))
    return out


def M_prime_bessel(y, smax=4000, dmax=60):
    """Direct (slow) Bessel series, for cross-checking M_prime at a few points."""
    y = np.atleast_1d(np.asarray(y, dtype=float))
    sy = np.sqrt(y)
    out = np.zeros_like(y)
    for d in range(1, dmax + 1):
        f = factor_sqfree(d)
        if f is None:
            continue
        w = _w(f)
        s = np.arange(1, smax + 1)
        s = s[np.gcd(s, d) == 1]
        out += w / d * (j1(4 * np.pi * np.outer(sy, s) / d) / s[None, :]).sum(axis=1)
    return (12 / np.pi) * sy * out


# ---------- kernel ----------
def G_kernel(y, t0, sig, k=2, nt=4001):
    t = np.linspace(t0 - 7 * sig, t0 + 7 * sig, nt)
    t = t[t > 0.05]
    H = np.exp(-0.5 * ((t - t0) / sig) ** 2)
    ph0 = 2 * np.imag(loggamma(k / 2 + 1j * t)) - 2 * t * np.log(2 * np.pi)
    dt = t[1] - t[0]
    y = np.atleast_1d(np.asarray(y, dtype=float))
    out = np.empty(y.shape, dtype=complex)
    for i0 in range(0, len(y), 2000):
        yy = y[i0:i0 + 2000]
        ph = ph0[None, :] - t[None, :] * np.log(yy)[:, None]
        out[i0:i0 + 2000] = (H[None, :] * np.exp(1j * ph)).sum(axis=1) * dt
    return out


def H_int(t0, sig):
    return sig * np.sqrt(2 * np.pi)


def H_hat(xi, t0, sig):
    return sig * np.sqrt(2 * np.pi) * np.exp(-0.5 * (sig * xi) ** 2) * np.exp(-1j * xi * t0)
