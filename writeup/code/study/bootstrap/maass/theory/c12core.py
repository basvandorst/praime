# Core of the weight-12 certificate (own code, python-flint arb; tau(n) from PARI's ramanujantau, checked against the
# q-expansion of Delta in validate_exact12.py).
#
# Lambda(s) = Gamma_C(s + 11/2) L(s, Delta) = 2 int_0^oo Delta(iy) y^{s+11/2} dy/y;  Xi(t) = Lambda(1/2+it) = int phi(u) e^{itu} du,
# phi(u) = 2 e^{6u} Delta(i e^u)  (even: Delta(i/y) = y^12 Delta(iy); positive: product formula).
# Closed form: (phi*phi)(w) = 8 e^{6w} sum_N (tau*tau)(N) K0(4 pi sqrt(N) e^{w/2}).
# With z = 4 pi sqrt(N) e^{w/2}, d/dw = (1/2) z d/dz, and D^j[z^12 K0] = z^12 (p_j K0 + q_j K1),
#   p_{j+1} = 12 p_j + z p_j' - z q_j,  q_{j+1} = 11 q_j + z q_j' - z p_j:
#   h_j(n) := (d/dw)^j (phi*phi)(log n) = 8 n^6 sum_N (tau*tau)(N) 2^{-j} [p_j(z) K0(z) + q_j(z) K1(z)],  z = 4 pi sqrt(N n).
# Tail N > Nmax: |(tau*tau)(N)| <= N^{11/2} d_4(N) <= 8 N^7 (Deligne |tau(n)| <= d(n) n^{11/2}; d_4 <= d^3 <= 8 N^{3/2}).
import sys; sys.dont_write_bytecode = True
import math
from fractions import Fraction
from flint import arb, ctx

_TAU = None
def tau_table(NT):
    global _TAU
    if _TAU is None or len(_TAU) <= NT:
        import cypari2
        pari = cypari2.Pari()
        _TAU = [0] + [int(pari.ramanujantau(n)) for n in range(1, NT + 1)]
    return _TAU
_TT = {}
def tautau(N):
    if N not in _TT:
        t = tau_table(max(N, 4000))
        _TT[N] = sum(t[d]*t[N//d] for d in range(1, N + 1) if N % d == 0)
    return _TT[N]

_PQ = {}
def pq(J, a=12):
    if (J, a) in _PQ: return _PQ[(J, a)]
    p, q = [[Fraction(1)]], [[Fraction(0)]]
    for j in range(J):
        A, B = p[-1], q[-1]; L = max(len(A), len(B)) + 1
        A = A + [Fraction(0)]*(L - len(A)); B = B + [Fraction(0)]*(L - len(B))
        na = [Fraction(0)]*L; nb = [Fraction(0)]*L
        for i in range(L):
            na[i] += a*A[i] + i*A[i]; nb[i] += (a - 1)*B[i] + i*B[i]
            if i + 1 < L: na[i + 1] -= B[i]; nb[i + 1] -= A[i]
        while len(na) > 1 and na[-1] == 0: na.pop()
        while len(nb) > 1 and nb[-1] == 0: nb.pop()
        p.append(na); q.append(nb)
    _PQ[(J, a)] = (p, q); return p, q
def farb(fr): return arb(fr.numerator)/fr.denominator

def _asym(z, nu, prec):
    mu = 4*nu*nu; term = arb(1); s = arb(0); k = 0; target = arb(2)**(-prec - 8)
    while True:
        s += term; k += 1
        nt = term*(mu - (2*k - 1)**2)/(8*k*z)
        if abs(nt).upper() < target: return (arb.pi()/(2*z)).sqrt()*(-z).exp()*(s + arb(0, abs(nt).upper()))
        if abs(nt).upper() > abs(term).upper() and k > 1: return None
        term = nt
_KC = {}
def K01(M):
    """K0, K1 at z = 4 pi sqrt(M), M a positive integer (asymptotic series with the DLMF 10.40(ii) bound, else arb)."""
    prec = ctx.prec
    if (M, prec) in _KC: return _KC[(M, prec)]
    z = 4*arb.pi()*arb(M).sqrt()
    a0 = _asym(z, 0, prec); a1 = _asym(z, 1, prec)
    if a0 is not None and a1 is not None:
        _KC[(M, prec)] = (a0, a1); return a0, a1
    extra = int(3*float(z.mid())) + 60
    while True:
        old = ctx.prec; ctx.prec = prec + extra
        try:
            zz = 4*arb.pi()*arb(M).sqrt(); k0 = zz.bessel_k(0); k1 = zz.bessel_k(1)
        finally:
            ctx.prec = old
        if all(k.rad() <= abs(k.mid())*arb(2)**(-prec + 4) for k in (k0, k1)):
            _KC[(M, prec)] = (+k0, +k1); return _KC[(M, prec)]
        extra *= 2; assert extra < 40000

def gamma_upper_bound(s, X):
    """Upper bound for Gamma(s, X), X > s - 1 > 0: X^{s-1} e^{-X} X/(X - s + 1)."""
    assert X > s - 1
    return X**(s - 1)*(-X).exp()*X/(X - s + 1)

def tail_N(n, J, PQbar, Nmax):
    """Rigorous bound for sum_{N > Nmax} 8 N^7 PQbar z^J sqrt(pi/(2z)) e^{-z + 1/(2z)}, z = 4 pi sqrt(N n)
    (z >= 1, so |p_j|,|q_j| <= PQbar z^J).  f(x) = A x^b e^{-c sqrt x} is decreasing for c sqrt x > 2b; then
    sum_{N > Nmax} f(N) <= int_{Nmax}^oo f = 2 A c^{-(2b+2)} Gamma(2b+2, c sqrt(Nmax))."""
    c = 4*arb.pi()*arb(n).sqrt()
    b = arb(7) + arb(J)/2 - arb(1)/4
    z0 = c*arb(Nmax).sqrt()
    assert z0 > 2*b + 2 and z0 > 2*b + 1
    A = 8*PQbar*c**J*(arb.pi()/(2*c)).sqrt()*(1/(2*z0)).exp()      # f(x) = A x^{7 + J/2 - 1/4} e^{-c sqrt x}
    return 2*A*c**(-(2*b + 2))*gamma_upper_bound(2*b + 2, z0)

def h_values(n, orders):
    """[h_j(n) for j in orders] (arb balls; the N-tail in the radius)."""
    J = max(orders); p, q = pq(J)
    pa = {j: [farb(x) for x in p[j]] for j in orders}; qa = {j: [farb(x) for x in q[j]] for j in orders}
    sums = {j: arb(0) for j in orders}
    Nmax = 1
    while True:
        z = 4*math.pi*math.sqrt(Nmax*n)
        if z > 2*(7 + J/2) + 20 and -z + J*math.log(z) + 8*math.log(Nmax + 1) < -(ctx.prec + 20)*math.log(2) - 10: break
        Nmax += 1
    for N in range(1, Nmax + 1):
        w = tautau(N)
        if w == 0: continue
        z = 4*arb.pi()*arb(N*n).sqrt(); k0, k1 = K01(N*n)
        pw = [arb(1)]
        for i in range(J + 1): pw.append(pw[-1]*z)
        for j in orders:
            pv = sum((pa[j][i]*pw[i] for i in range(len(pa[j]))), arb(0))
            qv = sum((qa[j][i]*pw[i] for i in range(len(qa[j]))), arb(0))
            sums[j] += w*(pv*k0 + qv*k1)/arb(2)**j
    pref = 8*arb(n)**6; out = []
    for j in orders:
        PQ = farb(sum(abs(x) for x in p[j]) + sum(abs(x) for x in q[j]))/arb(2)**j
        t = tail_N(n, j, PQ, Nmax)
        out.append(pref*sums[j] + arb(0, (pref*t).upper()))
    return out

def g_basis(n, D):
    hv = h_values(n, [2*k for k in range(D + 1)])
    return [(-1)**k*hv[k] for k in range(D + 1)]
