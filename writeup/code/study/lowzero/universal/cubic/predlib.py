"""predlib.py: predicted murmuration density of S3 cubic fields (family average of the one-point law with lines).

m_sgn(y) = E_K[ sqrt|D| a_K(p) ] at y = p/|D| (p unramified, fields with X <= |D| < 2X of sign sgn, X -> oo):
    m(y) = 2 pi sum_rho B(rho) K_sgn(4 pi sqrt(rho y)) - 1/(2 sqrt y),
    B(rho) = prod_l b_l(v_l(rho)),  sum_rho B(rho) rho^-s = zeta(3) zeta(2s) zeta(3s) prod_l h_l(l^-s),
    h_l(X) = 1 - (X^3 + X^4)/l + (X + X^5 - 1/X)/l^2.
Squares summed by Poisson: B = (squares) * w, w(a/b) = zeta(3) mu(b) b^-2 omega(a) [(a,b) = 1, b squarefree],
    omega multiplicative, omega(l^v) = coefficient of X^v in h_l(X)/(1 - X^3)  (v >= 0).
Bin integrals:  int_0^y m_total = 2 pi zeta(3) F(y),  F(y) = sum_{rho'} w(rho')/zeta(3) * C(rho' y)/rho'  (kerlib.C).
"""
import math, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kerlib as KL
from scipy.special import j1, y1, k1

ZETA3 = 1.2020569031595942854
ZETA2 = math.pi**2 / 6
ZP2_OVER_Z2 = -0.5699609930944961  # zeta'(2)/zeta(2)

def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for q in range(2, int(n**0.5) + 1):
        if s[q]: s[q * q::q] = False
    return np.nonzero(s)[0]

def h_coef(l):
    return {-1: -1.0 / l**2, 0: 1.0, 1: 1.0 / l**2, 3: -1.0 / l, 4: -1.0 / l, 5: 1.0 / l**2}

def omega_l(l, v):
    h = h_coef(l); return sum(h.get(v - 3 * j, 0.0) for j in range(0, v // 3 + 2))

def enumerate_a(Amax, tau):
    """All a <= Amax with omega(a) != 0 and |omega(a)| a^(-3/4) >= tau: list of (a, omega(a), prime set)."""
    P = primes_upto(int(Amax)); out = []
    def rec(i, a, w, ps):
        out.append((a, w, ps))
        for j in range(i, len(P)):
            l = int(P[j])
            if a * l > Amax or abs(w) * a**-0.75 * l**-2.25 < tau: break   # best case over v is v = 3 (|omega| <= 1 - 1/l)
            v, al = 1, a * l
            while al <= Amax:
                wl = omega_l(l, v)
                if wl != 0.0 and abs(w * wl) * al**-0.75 >= tau:      # |omega| a^(-3/4) only decreases along a branch
                    rec(j + 1, al, w * wl, ps + (l,))
                v += 1; al *= l
    rec(0, 1, 1.0, ())
    return [(a, w, ps) for (a, w, ps) in out if abs(w) * a**-0.75 >= tau]

_SQF = None
def squarefree_upto(n):
    global _SQF
    if _SQF is None or len(_SQF[0]) < n + 1:
        mu = np.ones(n + 1, np.int64); isq = np.ones(n + 1, bool)
        for q in primes_upto(n):
            mu[q::q] *= -1; isq[q * q::q * q] = False
        mu[~isq] = 0; mu[0] = 0
        _SQF = (mu, None)
    mu = _SQF[0][:n + 1]
    b = np.nonzero(mu)[0]
    return b, mu[b]

def F_of_y(edges, sgn, Amax=10**6, tau=1e-6, B0=300, Bcap=6000, verbose=False):
    """F(y) at each edge (array), sum over rho' = a/b; returns F and bookkeeping."""
    edges = np.asarray(edges, float); ymax = edges.max()
    alist = enumerate_a(Amax, tau)
    ball, muall = squarefree_upto(Bcap)
    F = np.zeros_like(edges); ntail = 0; nneg = 0
    for (a, w, ps) in alist:
        Bexp = int(min(Bcap, max(B0, 20 * a * ymax))) if a <= 100 else int(max(B0, min(Bcap, Bcap * 100.0 / a)))
        sel = ball <= Bexp; b = ball[sel].astype(float); mu = muall[sel].astype(float)
        if ps:
            cop = np.ones(len(b), bool)
            for l in ps: cop &= (ball[sel] % l != 0)
            b = b[cop]; mu = mu[cop]
        Z = a * edges[None, :] / b[:, None]
        Cv = np.zeros_like(Z); pos = Z > 0
        Cv[pos] = KL.C(Z[pos], sgn)
        F += w * np.sum((mu / b**2)[:, None] * (b[:, None] / a) * Cv, axis=0)
        if a * ymax / Bexp <= 0.05:   # exact tail b > Bexp from the small-Z expansion
            ntail += 1
            S32 = (1 / 2.6123753486854883) / np.prod([1 - l**-1.5 for l in ps]) if ps else 1 / 2.6123753486854883
            S2 = (1 / ZETA2) / np.prod([1 - l**-2.0 for l in ps]) if ps else 1 / ZETA2
            T32 = S32 - np.sum(mu * b**-1.5); T2 = S2 - np.sum(mu * b**-2.0)
            if sgn < 0:
                tail = np.sqrt(edges) / (2 * math.pi * math.sqrt(a)) * T32 - edges / 2 * T2
            else:
                sl = -ZP2_OVER_Z2 - sum(math.log(l) / (l * l - 1) for l in ps)
                S2ln = -S2 * sl
                T2ln = S2ln - np.sum(mu * b**-2.0 * np.log(b))
                with np.errstate(divide='ignore', invalid='ignore'):
                    lg = np.where(edges > 0, np.log(np.maximum(a * edges, 1e-300)), 0.0)
                tail = np.sqrt(edges) / (2 * math.pi * math.sqrt(a)) * T32 + edges / math.pi * ((lg - 1 + 2 * KL.EULER) * T2 - T2ln)
            F += w * tail
        else:
            nneg += 1
    if verbose: print('  F_of_y: %d values of a (Amax %g, tau %g); exact b-tail for %d, neglected for %d' % (len(alist), Amax, tau, ntail, nneg))
    return F

def predict_bins(edges, sgn, **kw):
    """Predicted bin averages of m_sgn(y) (primes only) and of the line sum without the prime-square term."""
    edges = np.asarray(edges, float); h = np.diff(edges)
    F = F_of_y(edges, sgn, **kw)
    tot = 2 * math.pi * ZETA3 * np.diff(F) / h
    sq = np.diff(np.sqrt(edges)) / h
    return tot - sq, tot

# ---- comparators -------------------------------------------------------------------------------------------
def K1prim(x, sgn):
    """antiderivative helper: int_0^y K(4 pi sqrt(rho y')) dy' = sqrt(y) K1prim(4 pi sqrt(rho y)) / (2 pi sqrt(rho))."""
    x = np.asarray(x, float); out = np.zeros_like(x); p = x > 0
    out[p] = j1(x[p]) if sgn < 0 else (-y1(x[p]) - 2 / math.pi * k1(x[p]))
    return out

def single_line_bins(edges, sgn, amp=ZETA3):
    """2 pi amp K(4 pi sqrt y), bin averages (the r = 1 line alone)."""
    edges = np.asarray(edges, float); h = np.diff(edges)
    G = 2 * math.pi * amp * np.sqrt(edges) * K1prim(4 * math.pi * np.sqrt(edges), sgn) / (2 * math.pi)
    return np.diff(G) / h

def gl2_prime_level_bins(edges, sgn, smax=4000, dmax=400):
    """GL(2) prime-level natural-weight lines rho = (s/d)^2, A = (6/pi^2) s^-2 prod_{p|d} p^2/((p-1)^2 (p+1)),
    transplanted to the kernel K_sgn (no nebentypus, no squares term): bin averages of 2 pi sum A sqrt(rho) K."""
    edges = np.asarray(edges, float); h = np.diff(edges)
    b, mu = squarefree_upto(dmax); G = np.zeros_like(edges)
    for d in b:
        d = int(d); fac = [q for q in primes_upto(d) if d % q == 0]
        g = np.prod([q * q / ((q - 1.0)**2 * (q + 1.0)) for q in fac]) if fac else 1.0
        s = np.arange(1, smax + 1);
        if fac:
            m = np.ones(len(s), bool)
            for q in fac: m &= (s % q != 0)
            s = s[m]
        s = s.astype(float); A = 6 / math.pi**2 * g / s**2; sr = s / d
        x = 4 * math.pi * sr[:, None] * np.sqrt(edges)[None, :]
        G += np.sum(A[:, None] * np.sqrt(edges)[None, :] * K1prim(x, sgn), axis=0)
    return np.diff(G) / h

def tt_bins(edges, sgn):
    """Taniguchi-Thorne / Roberts secondary term alone: c_sgn y^(-1/3)."""
    edges = np.asarray(edges, float); h = np.diff(edges)
    c = -1.02494734415798 if sgn < 0 else -1.77526087516441
    return c * 1.5 * np.diff(edges**(2.0 / 3)) / h
