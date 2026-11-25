"""Independent re-implementation (review of round 141). Own theta (Stirling series, checked against mpmath),
own Kronecker table (PARI), own window counts, own predictions, own regression and jackknife."""
import math, os
import numpy as np

REPO = '/home/bas/praime'
ZFILE = os.path.join(REPO, 'writeup/data/riemann/lzeros/r138_Fn.npz')
T1, T2, H = 100.0, 600.0, 0.02
W = 0.3 * np.arange(1, 17)
OUT = os.path.dirname(os.path.abspath(__file__))

# Bernoulli numbers B_2k for the Stirling series of ln Gamma
_B = [1 / 6, -1 / 30, 1 / 42, -1 / 30, 5 / 66, -691 / 2730, 7 / 6]


def lngamma_stirling(z):
    z = np.asarray(z, complex)
    s = (z - 0.5) * np.log(z) - z + 0.5 * math.log(2 * math.pi)
    zp = z.copy(); z2 = z * z
    for k, b in enumerate(_B, start=1):
        s = s + b / ((2 * k) * (2 * k - 1) * zp)
        zp = zp * z2
    return s


def theta(t, q, a):
    """theta_chi(t) = (t/2) ln(q/pi) + Im lnGamma((1/2 + a + it)/2), valid for t >= 20 (Stirling, 7 terms)."""
    t = np.asarray(t, float)
    return t / 2 * math.log(q / math.pi) + np.imag(lngamma_stirling(0.25 + a / 2 + 0.5j * t))


def grid(w):
    M = int(math.floor((T2 - w - T1) / H + 1e-9)) + 1
    return T1 + H * np.arange(M)


def dn(g, D, w, tg=None):
    q = abs(D); a = 0 if D > 0 else 1
    tg = grid(w) if tg is None else tg
    c = np.searchsorted(g, tg + w, 'right') - np.searchsorted(g, tg, 'right')
    return c - (theta(tg + w, q, a) - theta(tg, q, a)) / math.pi


def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]


def kron_pari(Ds, pmax):
    import cypari2
    pari = cypari2.Pari()
    pari.allocatemem(400 * 10 ** 6)
    P = primes_upto(pmax)
    X = np.zeros((len(Ds), len(P)), np.int8)
    for i, D in enumerate(Ds):
        v = pari('apply(p->kronecker(%d,p), primes([2,%d]))' % (D, pmax))
        X[i] = np.array([int(x) for x in v], np.int8)
    return P, X


def phi_ratio(k):
    r = 1.0; n = k; p = 2
    while p * p <= n:
        if n % p == 0:
            r *= p / (p - 1)
            while n % p == 0: n //= p
        p += 1
    if n > 1: r *= n / (n - 1)
    return r


def load_zeros():
    z = np.load(ZFILE)
    return {int(k[1:]): np.sort(z[k]) for k in z.files}
