"""Three-point function of the zeros of a primitive Dirichlet L-function of conductor k (bad-prime set S = primes | k),
from the ratios conjecture (Conrey-Snaith 2008, triple correlation of zeta, eq. for J*(a; b1, b2)), with the rule for
chi mod k: every prime p | k is deleted from every Euler product and prime sum, zeta(1 +- x) -> (k/phi(k)) L(1 +- x, chi_0),
zeta'/zeta -> L'/L(., chi_0), and l = log(t/2pi) -> L = log(k t/2pi).

Unfolded positions r (mean spacing 1); raw separations eps = 2 pi r / L.
  T3(r1, r2, r3) = R3 - 1 - sum_pairs (R2 - 1)  (connected, unit density)  = (1/L^3) * sum of the six J* terms
  J*(a; b1, b2) = Q(a+b1, a+b2) - P1(a+b1) P3(a, b1, b2) - P1(a+b2) P3(a, b2, b1)
  P1(x) = e^{-L x} Zc(x),  Zc(x) = (k/phi)^2 L(1+x,chi0) L(1-x,chi0) A_k(x)
  P3(a, b, c) = B1(a+b, a+c) + (L'/L)(1+a+c, chi0) - (L'/L)(1+c-b, chi0)
GUE limit: Zc -> -1/x^2, L'/L(1+x) -> -1/x, B1, Q -> 0; then sum/L^3 = 2 S12 S23 S31.
Factorial cumulant k[3](l) = int_{[0,l]^3} T3 = 6 int_0^l (l - s) F(s) ds, F(s) = int_0^s T3(0, t, s) dt = 2 int_0^{s/2}.
Third cumulant kappa3 = k[3] + 3 Sigma2 - 2 l.
"""
import os, sys, math
import numpy as np
from scipy.interpolate import CubicSpline
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', '..', 'bk'))
import bklib

TWO_PI = 2 * math.pi


def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]


class CachedGrid(bklib.Grid):
    def __init__(self):
        z = np.load(os.path.join(HERE, 'grid.npz'))
        self.eps = z['eps']; self.eps_max = float(self.eps[-1])
        self.raw = {k: z[k] for k in ('Z2', 'Z2m1', 'regd', 'Phi', 'Phim1', 'd2lnPhiDiag')}
        self.sp = {k: ((CubicSpline(self.eps, v.real), CubicSpline(self.eps, v.imag)) if np.iscomplexobj(v) else CubicSpline(self.eps, v))
                   for k, v in self.raw.items()}


G = CachedGrid()
_gz = np.load(os.path.join(HERE, 'gtab.npz')); GSP = (CubicSpline(_gz['eps'], _gz['g'].real), CubicSpline(_gz['eps'], _gz['g'].imag))
EULER = 0.5772156649015329


# ---------------------------------------------------------------- one-variable functions at x = i*eps (eps real, any sign)
def g_reg(e, S=()):
    """(L'/L)(1 + i e, chi0) + 1/(i e): regular part. g(-e) = conj g(e)."""
    a = np.abs(e); v = GSP[0](a) + 1j * GSP[1](a); v = np.where(e < 0, np.conj(v), v)
    for p in S:
        v = v + math.log(p) / (p ** (1 + 1j * e) - 1)
    return v


def Zc_reg(e, S=()):
    """Zc(i e) - 1/e^2 (regular part of (k/phi)^2 |L(1+ie,chi0)|^2 A_k(ie)), = conj(bklib rego(e)) for e > 0."""
    a = np.abs(e)
    Z2m1 = G('Z2m1', a); Phi = G('Phi', a); Phim1 = G('Phim1', a)
    Psi = np.ones_like(a, dtype=complex); Psim1 = np.zeros_like(a, dtype=complex)
    for p in S:
        fp = bklib.f_p(p, a); F = 1 + a ** 2 * fp
        Psim1 = Psim1 * F + fp; Psi = Psi * F
    rego = Z2m1 * Phi * Psi + Phim1 * Psi + Psim1
    return np.where(e > 0, np.conj(rego), rego)


# ---------------------------------------------------------------- prime sums of two variables
class PrimeSums:
    def __init__(self, PB=10 ** 5):
        self.P = primes_upto(PB).astype(float); self.lnP = np.log(self.P)

    def B1_all(self, e1, e2, chunk=200):
        """B1(i e1, i e2) summed over all p <= PB (Conrey-Snaith 2008; checked against the recipe in check_B1.py)."""
        out = np.empty(len(e1), complex); P = self.P; lp = self.lnP
        for j in range(0, len(e1), chunk):
            E1 = np.exp(-1j * np.outer(e1[j:j + chunk], lp)); E2 = np.exp(-1j * np.outer(e2[j:j + chunk], lp))   # p^{-x}, p^{-y}
            num = (1 - E1) * (1 - E1 - E2 + E2 / P) * lp
            den = (1 - E2 / (P * E1)) * (1 - E2 / P) * (1 - 2 / P + E1 / P) * P ** 2 * E1 / E2
            out[j:j + chunk] = (num / den).sum(1)
        return out

    def Q_all(self, e1, e2, chunk=200):
        """Q(i e1, i e2) = -sum_p ln^3 p / ((p^{1+x}-1)(p^{1+y}-1)), p <= PB."""
        out = np.empty(len(e1), complex); P = self.P; lp = self.lnP
        for j in range(0, len(e1), chunk):
            X = P * np.exp(1j * np.outer(e1[j:j + chunk], lp)); Y = P * np.exp(1j * np.outer(e2[j:j + chunk], lp))
            out[j:j + chunk] = (-(lp ** 3) / ((X - 1) * (Y - 1))).sum(1)
        return out


def b1_p(p, e1, e2):
    lp = math.log(p); E1 = np.exp(-1j * e1 * lp); E2 = np.exp(-1j * e2 * lp)
    return (1 - E1) * (1 - E1 - E2 + E2 / p) * lp / ((1 - E2 / (p * E1)) * (1 - E2 / p) * (1 - 2 / p + E1 / p) * p ** 2 * E1 / E2)


def q_p(p, e1, e2):
    lp = math.log(p); X = p * np.exp(1j * e1 * lp); Y = p * np.exp(1j * e2 * lp)
    return -(lp ** 3) / ((X - 1) * (Y - 1))


# ---------------------------------------------------------------- the six terms
TERMS = (('w', 'u', 'v'), ('v', 'w', 'u'), ('u', 'w', 'v'))   # J*(-i A; i B1, i B2); the other three are complex conjugates


def T3_parts(u, v, w, L, S=(), PS=None, B1all=None, Qall=None):
    """connected T3 at unfolded positions (arrays). Returns dict: gue, Q, P (P = arithmetic P-part minus its GUE limit),
    each already divided by L^3 and summed over the six terms (real). B1all/Qall: optional precomputed all-prime sums
    (lists of 6 arrays for B1 in the order T1(x1,x2), T1(x2,x1), T2..., and 3 arrays for Q)."""
    k = TWO_PI / L; pos = {'u': u * k, 'v': v * k, 'w': w * k}
    totQ = 0; totP = 0; totG = 0; ib = 0
    for it, (A, B1n, B2n) in enumerate(TERMS):
        al = -pos[A]; be1 = pos[B1n]; be2 = pos[B2n]
        e1 = al + be1; e2 = al + be2; ed = be2 - be1           # x1 = i e1, x2 = i e2, d = b2 - b1 = i ed
        # Q
        if Qall is not None: Qv = Qall[it]
        else: Qv = PS.Q_all(e1, e2)
        for p in S: Qv = Qv - q_p(p, e1, e2)
        # the two P-terms: (x1, x2, d) and (x2, x1, -d)
        for (f1, f2, fd) in ((e1, e2, ed), (e2, e1, -ed)):
            if B1all is not None: B = B1all[ib]
            else: B = PS.B1_all(f1, f2)
            ib += 1
            for p in S: B = B - b1_p(p, f1, f2)
            ph = np.exp(-1j * L * f1)                           # e^{-L x1}, x1 = i f1
            h = B + g_reg(f2, S) - g_reg(fd, S)                 # P3 - P3_GUE
            P3G = -1 / (1j * f2) + 1 / (1j * fd)
            totP = totP - ph * (Zc_reg(f1, S) * (P3G + h) + h / f1 ** 2)
            totG = totG - ph * (1 / f1 ** 2) * P3G
        totQ = totQ + Qv
    c = 2 / L ** 3
    return {'gue_terms': (c * totG).real, 'Q': (c * totQ).real, 'P': (c * totP).real}


def sinc(x): return np.sinc(x)


def T3_gue(u, v, w): return 2 * sinc(v - u) * sinc(w - v) * sinc(w - u)


# ---------------------------------------------------------------- quadrature for k[3](l)
class Quad:
    """s (= r3 - r1) in pieces of width WS with NS Gauss nodes; t in [0, s/2] in ceil((s/2)/WT) pieces of NT nodes."""

    def __init__(self, lmax=8.0, WS=0.25, NS=8, WT=0.5, NT=8):
        xs, ws = np.polynomial.legendre.leggauss(NS); xt, wt = np.polynomial.legendre.leggauss(NT)
        npc = int(round(lmax / WS)); self.npc = npc; self.WS = WS
        s_nodes = (np.arange(npc)[:, None] * WS + (WS / 2) * (xs[None, :] + 1)).ravel(); s_w = np.tile((WS / 2) * ws, npc)
        S_, T_, Wt_, I_ = [], [], [], []
        for i, (s, sw) in enumerate(zip(s_nodes, s_w)):
            half = s / 2; m = max(1, int(math.ceil(half / WT))); h = half / m
            t = (np.arange(m)[:, None] * h + (h / 2) * (xt[None, :] + 1)).ravel(); tw = np.tile((h / 2) * wt, m)
            S_.append(np.full(len(t), s)); T_.append(t); Wt_.append(2 * tw); I_.append(np.full(len(t), i))
        self.s_nodes, self.s_w = s_nodes, s_w
        self.S = np.concatenate(S_); self.T = np.concatenate(T_); self.Wt = np.concatenate(Wt_); self.I = np.concatenate(I_)
        self.lgrid = (np.arange(npc) + 1) * WS

    def kf3(self, vals):
        """vals: T3 at the nodes (u=0, v=T, w=S). Returns k[3](l) on lgrid."""
        F = np.bincount(self.I, weights=vals * self.Wt, minlength=len(self.s_nodes))      # F(s) at s nodes
        f = F.reshape(self.npc, -1); w = self.s_w.reshape(self.npc, -1); s = self.s_nodes.reshape(self.npc, -1)
        I0 = np.cumsum((f * w).sum(1)); I1 = np.cumsum((f * w * s).sum(1))
        return 6 * (self.lgrid * I0 - I1)


def kf3_Q_closed(l, L, primes, mmax=60):
    """the Q-part of k[3](l) in closed form: -(6/pi^3) sum_p sum_{m,n>=1} sin(m th) sin(n th) sin((m+n) th)/(m n (m+n) p^(m+n)),
    th = pi l ln p / L  (exact cube integral of the six Q terms). Summed over the given primes."""
    l = np.atleast_1d(np.asarray(l, float)); out = np.zeros(len(l))
    P = np.asarray(list(primes), float)
    if len(P) == 0: return out
    for s in range(2, 2 * mmax):                       # s = m + n; terms with p^-s < 1e-17 dropped
        Ps = P[P ** (-float(s)) >= 1e-17]
        if len(Ps) == 0: break
        th = math.pi * np.outer(np.log(Ps), l) / L; ps = Ps[:, None] ** float(s)
        for m in range(1, s):
            n = s - m
            out += (np.sin(m * th) * np.sin(n * th) * np.sin(s * th) / (m * n * s * ps)).sum(0)
    return -(6 / math.pi ** 3) * out


def sigma2_bk(L, S, lgrid):
    """Bogomolny-Keating number variance (bklib) at fixed L for bad set S, on lgrid (multiples of 0.025)."""
    full = bklib.sigma2_full(G, L, tuple(S)); idx = [int(round(l / bklib.W)) - 1 for l in lgrid]
    return full[idx]


def sigma2_gue(lgrid):
    v = bklib.sigma2_gue(); idx = [int(round(l / bklib.W)) - 1 for l in lgrid]
    return v[idx]
