"""pred4.py: murmuration density of S4 quartic fields from the family-averaged one-point law (Mellin-Barnes evaluation).

m(y) = E_K[sqrt|D| a_K(p)], y = p/|D|.  Line sum m_lines(y) = sum_rho B(rho) k(rho y), k(z) = z^(-1/2) K_3(z), Mellin transform
  int_0^oo y^(w-1) m_lines(y) dy = D(w) Kt(w),  Kt(w) = prod_j Gamma_R(w + mu_j)/Gamma_R(1 - w + mu_j)   (the gamma ratio of X(1/2+it)),
  D(w) = sum_rho B(rho) rho^-w = A1 zeta(2w) zeta(3w) zeta(4w) prod_l P_l(l^-w),   A1 = prod_l (1 + l^-2 - l^-3 - l^-4)^-1,
  P_l(X) = (1 - X^2/l) h_l(X), h_l(X) = 1 - X^-1/l^2 + X^2/l - X^4/l - X^5/l + X^5/l^2 - X^6/l + X^7/l^2.
Prime squares contribute (1/2) y^(-1/2) (E[sum alpha^2] = 1): m = m_lines - 1/(2 sqrt y).
Primitive M1(y) = int_0^y m = [res. at w = 1/2: y^(1/2), cancelled by the squares] + r3 y^(2/3)/(2/3) + (1/2 pi i) int_(c) D Kt y^(1-w)/(1-w) dw,
c = 0.3 (between the poles 1/4 and 1/3).  A generic Dirichlet series 'model' can replace D (comparators, cubic check, empirical local types).
"""
import math, numpy as np
from scipy.special import loggamma
import mpmath as mp

def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for q in range(2, int(n**0.5) + 1):
        if s[q]: s[q * q::q] = False
    return np.nonzero(s)[0]
PR = primes_upto(10**4)
BERN = [float(mp.bernoulli(2 * k)) / math.factorial(2 * k) for k in range(1, 15)]

def zeta_em(s, chunk=4000):
    """Riemann zeta for a complex array s (Re s >= 0.5), Euler-Maclaurin with N ~ |s|/pi + 30 and 14 Bernoulli terms."""
    s = np.asarray(s, complex); out = np.empty_like(s); idx = np.argsort(np.abs(s)); 
    for a in range(0, len(s), chunk):
        ii = idx[a:a + chunk]; ss = s[ii]; N = int(np.abs(ss).max() / math.pi) + 30
        n = np.arange(1, N, dtype=float); acc = np.zeros(len(ss), complex)
        for b in range(0, len(n), 512):
            nb = n[b:b + 512]; acc += np.exp(-np.outer(ss, np.log(nb))).sum(1)
        lnN = math.log(N); NS = np.exp(-ss * lnN)
        acc += NS / 2 + N * NS / (ss - 1)
        poch = ss.copy(); powN = NS / N    # (s)_{2k-1} N^{-s-2k+1}
        for k, B in enumerate(BERN):
            acc += B * poch * powN
            poch = poch * (ss + 2 * k + 1) * (ss + 2 * k + 2); powN = powN / N**2
        out[ii] = acc
    return out

def zeta_real(x):
    return float(mp.zeta(x))

LNPI = math.log(math.pi)
def Kt(w, mus):
    """gamma ratio prod_j Gamma_R(w + mu_j)/Gamma_R(1 - w + mu_j) for complex array w."""
    w = np.asarray(w, complex); lg = np.zeros_like(w)
    for m in mus:
        lg += -(w + m) / 2 * LNPI + loggamma((w + m) / 2) + (1 - w + m) / 2 * LNPI - loggamma((1 - w + m) / 2)
    return np.exp(lg)

SIG_MU = {0: (0, 0, 0), 1: (0, 0, 1), 2: (0, 1, 1)}     # r2 = number of complex places
# local types of etale quartic algebras over Q_l: name -> (mass as function of l, coefficients of 1/L_l(X))
TYPES = {'1111': (lambda l: 1 / 24, [1, -3, 3, -1]), '211': (lambda l: 1 / 4, [1, -1, -1, 1]), '22': (lambda l: 1 / 8, [1, 1, -1, -1]),
         '31': (lambda l: 1 / 3, [1, 0, 0, -1]), '4': (lambda l: 1 / 4, [1, 1, 1, 1]), '1^2 11': (lambda l: 1 / (2 * l), [1, -2, 1]),
         '1^2 2': (lambda l: 1 / (2 * l), [1, 0, -1]), '1^2 1^2': (lambda l: 1 / (2 * l**2), [1, -1]), '2^2': (lambda l: 1 / (2 * l**2), [1, 1]),
         '1^3 1': (lambda l: 1 / l**2, [1, -1]), '1^4': (lambda l: 1 / l**3, [1])}
def poly(c, X):
    return sum(ck * X**k for k, ck in enumerate(c))
def F_type(name, l, X):
    """per-field local generating function sum_v F(v) X^v = (l/(l-1)) L^-1(1/(lX)) L(X)."""
    c = TYPES[name][1]; return l / (l - 1) * poly(c, 1 / (l * X)) / poly(c, X)
def E_local(l, X, freqs=None):
    if freqs is None:
        Z = 1 + 1 / l + 2 / l**2 + 1 / l**3; return sum(TYPES[n][0](l) * F_type(n, l, X) for n in TYPES) / Z
    return sum(f * F_type(n, l, X) for n, f in freqs.items() if f)
def P_l(l, X):
    return (1 - X**2 / l) * (1 - 1 / (X * l**2) + X**2 / l - X**4 / l - X**5 / l + X**5 / l**2 - X**6 / l + X**7 / l**2)
A1 = float(np.prod(1 / (1 + PR**-2.0 - PR**-3.0 - PR**-4.0)) * math.exp(-sum(1.0 / q**2 for q in range(10**4 + 1, 10**6, 1) if all(q % p for p in (2, 3, 5, 7))) * 0))

def A1_exact():
    P = primes_upto(10**7).astype(float); return float(np.exp(-np.sum(np.log1p(P**-2 - P**-3 - P**-4))))

def Q_euler(w, L=10**4):
    """prod_l P_l(l^-w) = R(w)/zeta(2 - w), R = prod_{l<=L} P_l/(1 - l^(w-2))  (Re w = 0.3: tail of R ~ 1e-6)."""
    w = np.asarray(w, complex); lr = np.zeros_like(w)
    for l in PR[PR <= L]:
        X = np.exp(-w * math.log(l)); lr += np.log(P_l(float(l), X) / (1 - 1 / (X * l**2)))
    return np.exp(lr) / zeta_em(2 - w)

def Q_real(x, L=10**7):
    P = primes_upto(L).astype(float); X = P**-x
    return float(np.exp(np.sum(np.log(P_l(P, X)))))

class Model:
    """a Dirichlet series D(w) on the contour, with the residues to the right of the contour."""
    def __init__(self, c, t, Dvals, poles):   # poles: list of (w0, residue of D at w0)
        self.c, self.t, self.D, self.poles = c, t, Dvals, poles

def quartic_model(c=0.3, T=1500.0, dt=0.01, local_override=None):
    """D(w) of the S4 family. local_override: {l: callable X -> E_l^emp(X)} replaces the Bhargava local factor at l."""
    t = np.arange(0, T + dt / 2, dt); w = c + 1j * t
    Dv = A1_exact() * zeta_em(2 * w) * zeta_em(3 * w) * zeta_em(4 * w) * Q_euler(w)
    q3 = Q_real(1 / 3); res3 = A1_exact() * zeta_real(2 / 3) * zeta_real(4 / 3) * q3 / 3
    res2 = A1_exact() * zeta_real(1.5) * zeta_real(2.0) * Q_real(0.5) / 2
    if local_override:
        for l, Ef in local_override.items():
            Dv = Dv * Ef(np.exp(-w * math.log(l))) / E_local(l, np.exp(-w * math.log(l)))
            res3 *= Ef(l**(-1 / 3)) / E_local(l, l**(-1 / 3)); res2 *= Ef(l**-0.5) / E_local(l, l**-0.5)
    return Model(c, t, Dv, [(0.5, res2), (1 / 3, res3)])

def M1(model, mus, y, squares=True, window=None):
    """primitive int_0^y m for an array y: residues + contour integral (Simpson in t); squares=True subtracts the prime squares."""
    y = np.asarray(y, float); c, t = model.c, model.t; w = c + 1j * t
    g = model.D * Kt(w, mus) / (1 - w)
    if window is not None: g = g * window(t)
    wts = np.ones(len(t)); wts[1:-1:2] = 4; wts[2:-1:2] = 2; wts *= (t[1] - t[0]) / 3
    out = np.zeros(len(y))
    for a in range(0, len(y), 50):
        yy = y[a:a + 50]; ph = np.exp(np.outer(np.log(yy), 1 - w))           # y^(1-w)
        out[a:a + 50] = (ph * (g * wts)[None, :]).sum(1).real / math.pi
    for w0, r in model.poles:
        out += r * float(Kt(np.array([w0]), mus)[0].real) * y**(1 - w0) / (1 - w0)
    if squares: out -= np.sqrt(y)
    return out

def bin_means(model, mus, edges, **kw):
    e = np.asarray(edges, float); F = M1(model, mus, e, **kw); return np.diff(F) / np.diff(e)

def tail_coeffs(model, mus):
    """small-y tail of m: coefficients of y^(-1/2) (lines; cancelled by squares) and y^(-1/3)."""
    return {('%.4f' % w0): r * float(Kt(np.array([w0]), mus)[0].real) for w0, r in model.poles}

# ---- accelerated real Euler product, cubic-field model (check against the S3 table), single line -------------------------
def Q_real_acc(x, L=10**6):
    """prod_l P_l(l^-x) = prod_l [P_l/(1 - l^(x-2))] / zeta(2 - x)."""
    Pp = primes_upto(L).astype(float); X = Pp**-x
    return float(np.exp(np.sum(np.log(P_l(Pp, X) / (1 - Pp**(x - 2))))) / zeta_real(2 - x))

def h_cub(l, X):
    return 1 - (X**3 + X**4) / l + (X + X**5 - 1 / X) / l**2

def cubic_model(c=0.3, T=1500.0, dt=0.01):
    """D(w) of the S3 cubic-field family (round 161): zeta(3) zeta(2w) zeta(3w) prod_l h_l(l^-w); acceleration 1/zeta(2-w), 1/zeta(1+3w)."""
    t = np.arange(0, T + dt / 2, dt); w = c + 1j * t; lr = np.zeros_like(w)
    for l in PR:
        X = np.exp(-w * math.log(l)); lr += np.log(h_cub(float(l), X) / ((1 - 1 / (X * l**2)) * (1 - X**3 / l)))
    z3 = zeta_real(3.0)
    Dv = z3 * zeta_em(2 * w) * zeta_em(3 * w) * np.exp(lr) / (zeta_em(2 - w) * zeta_em(1 + 3 * w))
    Pp = primes_upto(10**6).astype(float)
    def qr(x): return float(np.exp(np.sum(np.log(h_cub(Pp, Pp**-x) / ((1 - Pp**(x - 2)) * (1 - Pp**(-1 - 3 * x)))))) / (zeta_real(2 - x) * zeta_real(1 + 3 * x)))
    return Model(c, t, Dv, [(0.5, z3 * zeta_real(1.5) * qr(0.5) / 2), (1 / 3, z3 * zeta_real(2 / 3) * qr(1 / 3) / 3)])

def single_line_model(amp, c=0.3, T=1500.0, dt=0.01):
    t = np.arange(0, T + dt / 2, dt); return Model(c, t, amp * np.ones(len(t), complex), [])

CUBIC_MU = {-1: (0, 1), 1: (0, 0)}     # Gamma_C (D < 0), Gamma_R^2 (D > 0)
