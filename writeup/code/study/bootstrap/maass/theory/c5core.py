# Core of the chi_5 certificate (even character; own code, python-flint arb).
#
# Lambda(s) = 5^{s/2} Gamma_R(s) L(s, chi_5),  chi_5 = (n/5);  theta0(x) = sum_{n>=1} chi(n) exp(-pi n^2 x/5),
# theta0(1/x) = x^{1/2} theta0(x) (root number +1);  Lambda(s) = int_0^oo theta0(x) x^{s/2} dx/x;
# Xi(t) = Lambda(1/2+it) = int phi(u) e^{itu} du,  phi(u) = 2 e^{u/2} theta0(e^{2u})  (even).
# Closed form: (phi*phi)(w) = 4 e^{w/2} sum_{N>=1} chi(N) d(N) K0(c N e^w),  c = 2 pi/5.
# D = z d/dz:  D^j[z^{1/2} K0] = z^{1/2}(p_j K0 + q_j K1),  p_{j+1} = (1/2)p_j + z p_j' - z q_j,  q_{j+1} = -(1/2)q_j + z q_j' - z p_j.
# h_j(n) = (d/dw)^j (phi*phi)(log n) = 4 n^{1/2} sum_N chi(N) d(N) [p_j(z) K0(z) + q_j(z) K1(z)],  z = c N n.
import sys; sys.dont_write_bytecode = True
import math
from fractions import Fraction
from flint import arb, ctx

def chi(n): return (0, 1, -1, -1, 1)[n % 5]
def ndiv(n):
    c, d = 0, 1
    while d*d <= n:
        if n % d == 0: c += 1 if d*d == n else 2
        d += 1
    return c
_PQ = {}
def pq(J):
    if J in _PQ: return _PQ[J]
    p, q = [[Fraction(1)]], [[Fraction(0)]]
    for j in range(J):
        a, b = p[-1], q[-1]; L = max(len(a), len(b)) + 1
        a = a + [Fraction(0)]*(L - len(a)); b = b + [Fraction(0)]*(L - len(b))
        na = [Fraction(0)]*L; nb = [Fraction(0)]*L
        for i in range(L):
            na[i] += Fraction(1, 2)*a[i] + i*a[i]; nb[i] += Fraction(-1, 2)*b[i] + i*b[i]
            if i + 1 < L: na[i + 1] -= b[i]; nb[i + 1] -= a[i]
        while len(na) > 1 and na[-1] == 0: na.pop()
        while len(nb) > 1 and nb[-1] == 0: nb.pop()
        p.append(na); q.append(nb)
    _PQ[J] = (p, q); return p, q
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
    """K0, K1 at z = 2 pi M/5: asymptotic series (DLMF 10.40(ii) bound) or arb's bessel_k at raised precision."""
    prec = ctx.prec
    if (M, prec) in _KC: return _KC[(M, prec)]
    z = 2*arb.pi()*M/5
    a0 = _asym(z, 0, prec); a1 = _asym(z, 1, prec)
    if a0 is not None and a1 is not None:
        _KC[(M, prec)] = (a0, a1); return a0, a1
    extra = int(3*1.2567*M) + 60
    while True:
        old = ctx.prec; ctx.prec = prec + extra
        try:
            zz = 2*arb.pi()*M/5; k0 = zz.bessel_k(0); k1 = zz.bessel_k(1)
        finally:
            ctx.prec = old
        if all(k.rad() <= abs(k.mid())*arb(2)**(-prec + 4) for k in (k0, k1)):
            _KC[(M, prec)] = (+k0, +k1); return _KC[(M, prec)]
        extra *= 2; assert extra < 40000

def h_values(n, orders):
    """[h_j(n) for j in orders], balls, the N-tail in the radius (|p_j|,|q_j| <= sum|coef| z^j, K_nu <= sqrt(pi/2z) e^{-z+1/(2z)}, d(N) <= N)."""
    J = max(orders); p, q = pq(J); c = 2*arb.pi()/5
    pa = {j: [farb(x) for x in p[j]] for j in orders}; qa = {j: [farb(x) for x in q[j]] for j in orders}
    sums = {j: arb(0) for j in orders}
    Nmax = 1
    while True:
        z = 2*math.pi/5*Nmax*n
        if z > 2*J + 10 and -z + J*math.log(z) + 3*math.log(Nmax + 1) < -(ctx.prec + 20)*math.log(2) - 10: break
        Nmax += 1
    for N in range(1, Nmax + 1):
        ch = chi(N)
        if ch == 0: continue
        z = c*N*n; k0, k1 = K01(N*n)
        pw = [arb(1)]
        for i in range(J + 1): pw.append(pw[-1]*z)
        w = ch*ndiv(N)
        for j in orders:
            pv = sum((pa[j][i]*pw[i] for i in range(len(pa[j]))), arb(0))
            qv = sum((qa[j][i]*pw[i] for i in range(len(qa[j]))), arb(0))
            sums[j] += w*(pv*k0 + qv*k1)
    pref = 4*arb(n).sqrt(); out = []
    for j in orders:
        PQ = farb(sum(abs(x) for x in p[j]) + sum(abs(x) for x in q[j]))
        N = Nmax + 1; zN = c*N*n
        t = arb(N)*PQ*zN**j*(arb.pi()/(2*zN)).sqrt()*(-zN + 1/(2*zN)).exp()
        rho = (arb(N + 1)/N)**(j + 2)*(-c*n).exp(); assert rho.upper() < 1
        out.append(pref*sums[j] + arb(0, (pref*t/(1 - rho)).upper()))
    return out

def g_basis(n, D):
    hv = h_values(n, [2*k for k in range(D + 1)])
    return [(-1)**k*hv[k] for k in range(D + 1)]
