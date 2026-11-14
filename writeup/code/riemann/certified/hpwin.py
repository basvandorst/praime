"""Rigorous high-precision samples for the certified decoder, in ball arithmetic (Arb, via python-flint).

Nothing here is trusted beyond Arb and the two analytic bounds of cert.py:
- the zeros come from acb.zeta_zeros (Arb isolates them and checks the count with Turing's method), as balls;
- the kernel and its transform are entire functions written with 0F1, so no square root or branch is involved:
    K(v)    = beta^m (1-u^2)^m F_m(beta^2 (1-u^2)),        u = v/a, |v| < a
    Khat(w) = a sqrt(2 pi) beta^m F_nu(beta^2 - (a w)^2), nu = m + 1/2
    F_s(z)  = I_s(sqrt z)/sqrt(z)^s = 0F1(; s+1; z/4) / (2^s Gamma(s+1));
- the trivial zeros are integrated rigorously (acb.integral; the integrand is meromorphic, no poles on the path);
- the zeros above T are bounded by cert.tail_bound and cert.offline_bound (double precision, inflated by 1e-6 relative).
sample(x, a, beta, m, T) returns (zero side as a ball, delta as an upper bound).
"""
import math
from flint import arb, acb, ctx
import numpy as np
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cert import tail_bound, offline_bound

ctx.prec = int(os.environ.get('HPPREC', '192'))
_Z = []

def zeros_upto(T):
    """Imaginary parts of the zeros with 0 < gamma <= T, as balls (rigorous)."""
    while not _Z or _Z[-1].mid() <= T + 1:
        n0 = len(_Z) + 1
        for z in acb.zeta_zeros(n0, 200):
            _Z.append(z.imag)
    out = []
    for g in _Z:
        if g.upper() <= T: out.append(g)
        elif g.lower() > T: break
        else: raise ValueError(f'zero {g} straddles T = {T}')
    return out

def F(s, z):
    """I_s(sqrt z)/sqrt(z)^s, entire in z."""
    s = arb(s)
    return (z / 4).hypgeom_0f1(s + 1) / (arb(2) ** s * (s + 1).gamma())

def Fc(s, z):
    s = acb(s)
    return (z / 4).hypgeom_0f1(s + 1) / (acb(2) ** s * (s + 1).gamma())

def K(v, a, beta, m):
    a = arb(a); beta = arb(beta); v = arb(v)
    w = 1 - (v / a) ** 2
    if w.upper() <= 0: return arb(0)
    if w.lower() <= 0: raise ValueError('v at the edge of the support')
    return beta ** m * w ** m * F(m, beta * beta * w)

def K0(beta, m):
    beta = arb(beta); return beta ** m * F(m, beta * beta)

def Khat(w, a, beta, m):
    a = arb(a); beta = arb(beta)
    return a * (2 * arb.pi()).sqrt() * beta ** m * F(m + arb(1) / 2, beta * beta - (a * w) ** 2)

def P(x, a, beta, m):
    a = arb(a); beta = arb(beta); x = arb(x)
    return a * (2 * arb.pi()).sqrt() * beta ** m * F(m + arb(1) / 2, beta * beta + (a / (2 * x)) ** 2)

def triv(x, a, beta, m):
    """sqrt(x) int_{-a}^{a} K(v) sqrt(t) / (x (t^3 - t)) dv, t = x e^{v/x}; rigorous."""
    xa = acb(x); aa = acb(a); bb = acb(beta)
    def f(v, analytic):
        w = 1 - (v / aa) ** 2
        k = bb ** m * w ** m * Fc(m, bb * bb * w)
        t = xa * (v / xa).exp()
        return k * t.sqrt() / (xa * (t ** 3 - t))
    I = acb.integral(f, -a, a)
    return (arb(x).sqrt() * I.real)

def zero_side(x, a, beta, m, T, zs=None):
    zs = zeros_upto(T) if zs is None else zs
    xa = arb(x); L = xa.log(); s = arb(0)
    for g in zs:
        s += (g * L).cos() * Khat(g / xa, a, beta, m)
    return P(x, a, beta, m) - 2 / xa.sqrt() * s - triv(x, a, beta, m)

def delta(x, a, beta, m, T, NT):
    """Upper bound on |zero side computed from the zeros up to T - prime side| beyond the ball radius:
    zeros above T (Rosser's bound on N(t), up to 3e12) and above 3e12 (integration by parts)."""
    d = (2 / math.sqrt(x)) * tail_bound(x, T, NT, a, beta, m) + offline_bound(x, a, beta, m)
    return d * (1 + 1e-6) + 1e-300

def sample(x, a, beta, m, T):
    zs = zeros_upto(T)
    return zero_side(x, a, beta, m, T, zs), delta(x, a, beta, m, T, len(zs))
