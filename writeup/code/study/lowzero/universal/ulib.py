"""universal one-point law: shared library.
Analytic normalisation: L(s) = sum lambda(n) n^-s, Lambda(s) = Q^{s/2} prod_j Gamma_R(s + mu_j) L(s) = eps * conj-Lambda(1 - s).
theta_L(t) = (t/2) ln(Q/pi^d) + sum_j Im lnGamma((1/2 + mu_j + i t)/2), so X(1/2+it) = e^{-2 i theta}, L(1/2+it) = eps e^{-2i theta} conj L(1/2+it).
PARI (arithmetic) data [N, k, Vga]: mu_j = Vga_j + (k-1)/2, analytic s = arithmetic s - (k-1)/2."""
import math, cmath
import numpy as np
from scipy.special import loggamma

def theta(t, Q, mu):
    t = np.asarray(t, float)
    out = t / 2 * math.log(Q / math.pi ** len(mu))
    for m in mu:
        out = out + np.imag(loggamma((0.5 + m + 1j * t) / 2))
    return out

def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]

# ---------------------------------------------------------------- the L-functions of this study
# key: (builder, args, T2 window end, family)
CATALOGUE = [
    # GL(2) holomorphic, trivial character, level 1, weights 12-26 (rational, except weight 24: two embeddings)
    ('D12', ('mf', 1, 12, 1, 0, None), 600, 'GL2 level 1'),
    ('F16', ('mf', 1, 16, 1, 0, None), 600, 'GL2 level 1'),
    ('F18', ('mf', 1, 18, 1, 0, None), 600, 'GL2 level 1'),
    ('F20', ('mf', 1, 20, 1, 0, None), 600, 'GL2 level 1'),
    ('F22', ('mf', 1, 22, 1, 0, None), 600, 'GL2 level 1'),
    ('F26', ('mf', 1, 26, 1, 0, None), 600, 'GL2 level 1'),
    ('F24a', ('mf', 1, 24, 1, 0, 0), 600, 'GL2 level 1'),
    ('F24b', ('mf', 1, 24, 1, 0, 1), 600, 'GL2 level 1'),
    # GL(2) weight 2 with non-rational (real quadratic) coefficients: both embeddings
    ('W23a', ('mf', 23, 2, 1, 0, 0), 400, 'GL2 weight 2 irrational'),
    ('W23b', ('mf', 23, 2, 1, 0, 1), 400, 'GL2 weight 2 irrational'),
    ('W29a', ('mf', 29, 2, 1, 0, 0), 400, 'GL2 weight 2 irrational'),
    ('W29b', ('mf', 29, 2, 1, 0, 1), 400, 'GL2 weight 2 irrational'),
    ('W31a', ('mf', 31, 2, 1, 0, 0), 400, 'GL2 weight 2 irrational'),
    ('W31b', ('mf', 31, 2, 1, 0, 1), 400, 'GL2 weight 2 irrational'),
    # non-self-dual: weight 2, level 13, character of order 6 (f and its dual f-bar)
    ('N13', ('mf', 13, 2, 4, 0, None), 400, 'GL2 nebentypus (non-self-dual)'),
    ('N13d', ('mf', 13, 2, 10, 0, None), 400, 'GL2 nebentypus (non-self-dual)'),
    # degree 1 non-self-dual controls (theorem): chi mod 7 of order 3, chi mod 5 of order 4, and conjugates
    ('X7', ('dir', 7, 2), 1000, 'GL1 complex character'),
    ('X7d', ('dir', 7, 4), 1000, 'GL1 complex character'),
    ('X5', ('dir', 5, 2), 1000, 'GL1 complex character'),
    ('X5d', ('dir', 5, 3), 1000, 'GL1 complex character'),
    # GL(3): symmetric squares of elliptic curves
    ('S11a', ('sym', [0, -1, 1, -10, -20], 2), 200, 'GL3 Sym2 E'),
    ('S14a', ('sym', [1, 0, 1, 4, -6], 2), 200, 'GL3 Sym2 E'),
    ('S15a', ('sym', [1, 1, 1, -10, -10], 2), 200, 'GL3 Sym2 E'),
    ('S17a', ('sym', [1, -1, 1, -1, -14], 2), 200, 'GL3 Sym2 E'),
    ('S19a', ('sym', [0, 1, 1, -9, -15], 2), 200, 'GL3 Sym2 E'),
    ('S37a', ('sym', [0, 0, 1, -1, 0], 2), 200, 'GL3 Sym2 E'),
    # Artin (weight 1, S3) and the Dedekind zeta functions they sit in (products: zeta * L(rho))
    ('R23', ('artin', 'x^3-x-1'), 500, 'GL2 Artin S3'),
    ('R31', ('artin', 'x^3+x-1'), 500, 'GL2 Artin S3'),
    ('Z23', ('nf', 'x^3-x-1'), 300, 'product zeta*L(rho)'),
    ('Z31', ('nf', 'x^3+x-1'), 300, 'product zeta*L(rho)'),
    # GL(4): genus-2 curves (End = Z expected) and X1(13) (Jacobian of GL2 type: L(f) L(f-bar), a product)
    ('G249', ('g2', 'x^2+x', 'x^3+1'), 200, 'GL4 genus 2'),
    ('G277', ('g2', 'x^5-2*x^3+2*x^2-x', '1'), 200, 'GL4 genus 2'),
    ('J13', ('g2', 'x^5+x^4', 'x^3+x+1'), 200, 'product L(f)L(fbar) genus 2'),
]

def build(pari, spec):
    kind = spec[0]
    if kind == 'mf':
        _, N, k, chi, idx, emb = spec
        mf = pari('mfinit([%d,%d,Mod(%d,%d)],0)' % (N, k, chi, max(N, 1))) if chi != 1 else pari('mfinit([%d,%d],0)' % (N, k))
        B = pari.mfeigenbasis(mf)
        Lv = pari.lfunmf(mf, B[idx])
        if emb is not None:
            return Lv[emb]
        return Lv
    if kind == 'dir':
        _, q, m = spec
        return pari('lfuncreate([G=znstar(%d,1), znconreychar(G,%d)])' % (q, m))
    if kind == 'sym':
        _, a, m = spec
        return pari.lfunsympow(pari.ellinit(a), m)
    if kind == 'artin':
        return pari('lfundiv(lfuncreate(%s), lfuncreate(1))' % spec[1])
    if kind == 'nf':
        return pari('lfuncreate(%s)' % spec[1])
    if kind == 'g2':
        return pari('lfungenus2([%s,%s])' % (spec[1], spec[2]))
    raise ValueError(kind)

def params(pari, L):
    N, k, V = pari.lfunparams(L)
    k = float(k); mu = [float(v) + (k - 1) / 2 for v in V]
    return float(N), k, mu

def to_c(x):
    try:
        return complex(x)
    except Exception:
        return complex(x.lift()) if hasattr(x, 'lift') else complex(float(x))

def local_poly(pari, L, p, k):
    """analytic local polynomial coefficients [1, c1, ..., cd] with L_p = 1/P(p^-s)."""
    e = pari.lfuneuler(L, p)
    P = 1 / e
    co = [to_c(c) for c in pari.Vecrev(P)] if str(pari.type(P)) == 't_POL' else [to_c(P)]
    c0 = co[0]
    co = [c / c0 for c in co]
    return [c * p ** (-j * (k - 1) / 2) for j, c in enumerate(co)]

def series_inv(P, J):
    """coefficients of 1/P(X) up to X^J."""
    lam = [0j] * (J + 1); lam[0] = 1 / P[0]
    for j in range(1, J + 1):
        s = 0j
        for i in range(1, min(j, len(P) - 1) + 1):
            s += P[i] * lam[j - i]
        lam[j] = -s / P[0]
    return lam

def local_data(P, p, J=80):
    """F_p(1/p), S_p(v) for v = -2..2 (S_p(v) = sum_{m-h=v} conj(lambda(p^m)) mu(p^h) p^{-(m+h)/2})."""
    lam = series_inv(P, J)
    mu = P
    X = 1.0 / p
    F = sum(np.conj(lam[j]) * mu[j] * X ** j for j in range(len(mu)))
    S = {}
    for v in (-2, -1, 1, 2):
        s = 0j
        for h in range(len(mu)):
            m = h + v
            if 0 <= m <= J:
                s += np.conj(lam[m]) * mu[h] * p ** (-(m + h) / 2)
        S[v] = s
    return complex(F), S, lam


class Local:
    """analytic local polynomials P_p (L_p = 1/P_p(p^-s)) for every p in a list, and lambda(p) for all p <= Pmax."""
    def __init__(self, pari, L, spec, Q, k, Pmax):
        self.pari, self.L, self.spec, self.Q, self.k = pari, L, spec, Q, k
        an = pari.lfunan(L, Pmax)
        self.lam = np.array([to_c(an[n - 1]) for n in range(1, Pmax + 1)]) / np.arange(1, Pmax + 1) ** ((k - 1) / 2)
        self.Pmax = Pmax
        self.chi2 = None
        if spec[0] == 'mf' and spec[3] != 1:
            # character value at 2 (a primitive root mod 13 here): chi(2) = lambda(2)^2 - lambda(4)
            assert spec[1] == 13
            self.chi2 = self.lam[1] ** 2 - self.lam[3]
    def chi(self, p):
        if self.spec[0] != 'mf': raise ValueError
        N = self.spec[1]
        if N % p == 0: return 0
        if self.spec[3] == 1: return 1
        e = int(self.pari.znlog(p, self.pari.Mod(2, N)))
        return self.chi2 ** e
    def poly(self, p):
        if self.spec[0] == 'mf':
            N = self.spec[1]; lp = self.lam[p - 1]
            if N % p == 0: return [1, -lp]
            return [1, -lp, self.chi(p)]
        return local_poly(self.pari, self.L, p, self.k)
