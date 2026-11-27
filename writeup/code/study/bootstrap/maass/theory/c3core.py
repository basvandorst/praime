# Core of the chi_{-3} certificate (own code, python-flint arb ball arithmetic).
#
# Lambda(s) = 3^{s/2} Gamma_R(s+1) L(s, chi_{-3}),  Xi(t) = Lambda(1/2 + it) = int phi(u) e^{itu} du,
# phi(u) = (2/sqrt 3) e^{3u/2} sum_n n chi(n) exp(-pi n^2 e^{2u}/3)   (even in u: theta transformation).
# Closed form (derived in exact3.md, section 2):
#   (phi*phi)(w) = (4/3) e^{3w/2} sum_{N>=1} chi(N) d(N) N K0(c N e^w),   c = 2 pi / 3,  every real w.
# With D = z d/dz and z = c N e^w:  D^j [z^{3/2} K0(z)] = z^{3/2} [p_j(z) K0(z) + q_j(z) K1(z)],
#   p_{j+1} = (3/2) p_j + z p_j' - z q_j,   q_{j+1} = (1/2) q_j + z q_j' - z p_j,   p_0 = 1, q_0 = 0.
# So at w = log n (n a positive integer, e^w = n exactly):
#   h_j(n) := (d/dw)^j (phi*phi)(log n) = (4/3) n^{3/2} sum_N chi(N) d(N) N [p_j(z) K0(z) + q_j(z) K1(z)],  z = c N n.
# The function with transform Xi(t)^2 t^{2k} is g_k = (-1)^k h_{2k}.
from fractions import Fraction
from flint import arb, ctx

def chi(n):
    r = n % 3
    return 0 if r == 0 else (1 if r == 1 else -1)

def ndiv(n):
    c, d = 0, 1
    while d*d <= n:
        if n % d == 0:
            c += 1 if d*d == n else 2
        d += 1
    return c

_PQ = {}
def pq(J):
    """Exact coefficient lists (Fractions, index = power of z) of p_j, q_j for j <= J."""
    if J in _PQ: return _PQ[J]
    p, q = [[Fraction(1)]], [[Fraction(0)]]
    for j in range(J):
        a, b = p[-1], q[-1]; L = max(len(a), len(b)) + 1
        a = a + [Fraction(0)]*(L - len(a)); b = b + [Fraction(0)]*(L - len(b))
        na = [Fraction(0)]*L; nb = [Fraction(0)]*L
        for i in range(L):
            na[i] += Fraction(3, 2)*a[i] + i*a[i]; nb[i] += Fraction(1, 2)*b[i] + i*b[i]
            if i + 1 < L:
                na[i + 1] -= b[i]; nb[i + 1] -= a[i]
        while len(na) > 1 and na[-1] == 0: na.pop()
        while len(nb) > 1 and nb[-1] == 0: nb.pop()
        p.append(na); q.append(nb)
    _PQ[J] = (p, q)
    return p, q

def farb(fr):
    return arb(fr.numerator)/fr.denominator   # denominators are powers of 2: exact

def _asym(z, nu, prec):
    """K_nu(z) for real z > 0 from the asymptotic series with the DLMF 10.40(ii) bound (nu real, z > 0, l >= nu - 1/2:
    the remainder after l terms is at most the first neglected term in size).  Returns None if not accurate enough."""
    mu = 4*nu*nu; term = arb(1); s = arb(0); k = 0; target = arb(2)**(-prec - 8)
    while True:
        s += term
        k += 1
        nt = term*(mu - (2*k - 1)**2)/(8*k*z)
        if abs(nt).upper() < target and k >= nu:
            return (arb.pi()/(2*z)).sqrt()*(-z).exp()*(s + arb(0, abs(nt).upper()))
        if abs(nt).upper() > abs(term).upper() and k > 1:
            return None
        term = nt

def K01(M):
    """K0(z), K1(z) at z = 2 pi M / 3 (M a positive integer): asymptotic series when it reaches the working precision,
    otherwise arb's bessel_k at a raised precision (arb loses about 2 z log2(e) bits there; z is recomputed at the raised
    precision so that its own radius does not spread), repeated until the relative radius is small."""
    prec = ctx.prec
    z = 2*arb.pi()*M/3
    a0 = _asym(z, 0, prec); a1 = _asym(z, 1, prec)
    if a0 is not None and a1 is not None:
        return a0, a1
    zf = float(z.mid()); extra = int(3*zf) + 60
    while True:
        old = ctx.prec; ctx.prec = prec + extra
        try:
            zz = 2*arb.pi()*M/3
            k0 = zz.bessel_k(0); k1 = zz.bessel_k(1)
        finally:
            ctx.prec = old
        ok = all(float(k.rad()) <= float(abs(k.mid()))*2.0**(-prec + 4) for k in (k0, k1))
        if ok:
            return +k0, +k1
        extra *= 2
        assert extra < 20000, 'bessel_k failed'

def coef_abs_sum(cl):
    return sum(abs(c) for c in cl)

def h_values(n, orders, Nmax=None, tol_bits=None):
    """[h_j(n) for j in orders] as arb balls, with a rigorous bound for the tail N > Nmax added to the radius.
    Tail bound: |p_j(z)| <= Pbar_j z^j, |q_j(z)| <= Qbar_j z^j for z >= 1 (deg p_j, q_j <= j), K_nu(z) <= sqrt(pi/(2z)) e^{-z+1/(2z)}
    (from K_nu(z) = int_0^oo e^{-z cosh s} cosh(nu s) ds and cosh s >= 1 + s^2/2), d(N) <= N."""
    J = max(orders); p, q = pq(J)
    c = 2*arb.pi()/3
    if tol_bits is None: tol_bits = ctx.prec + 20
    pa = {j: [farb(x) for x in p[j]] for j in orders}; qa = {j: [farb(x) for x in q[j]] for j in orders}
    sums = {j: arb(0) for j in orders}
    N = 0
    if Nmax is None:
        # choose Nmax: z = cNn beyond max(2J, ...) and e^{-z} z^J small
        Nmax = 1
        while True:
            z = float(c.mid())*Nmax*n
            import math
            lg = -z + J*math.log(max(z, 1.0)) + 3*math.log(Nmax + 1)
            if z > 2*J + 10 and lg < -tol_bits*math.log(2) - 10: break
            Nmax += 1
    for N in range(1, Nmax + 1):
        ch = chi(N)
        if ch == 0: continue
        z = c*N*n; k0, k1 = K01(N*n)
        pw = [arb(1)]
        for i in range(J + 1): pw.append(pw[-1]*z)
        w = ch*ndiv(N)*N
        for j in orders:
            pv = sum((pa[j][i]*pw[i] for i in range(len(pa[j]))), arb(0))
            qv = sum((qa[j][i]*pw[i] for i in range(len(qa[j]))), arb(0))
            sums[j] += w*(pv*k0 + qv*k1)
    # rigorous tail N > Nmax
    import math
    out = []
    pref = arb(4)/3*arb(n)**arb(1.5) if False else arb(4)/3*arb(n)*arb(n).sqrt()
    for j in orders:
        Pb = sum(abs(x) for x in p[j]); Qb = sum(abs(x) for x in q[j])
        PQ = farb(Fraction(Pb) + Fraction(Qb))
        tail = arb(0); N = Nmax + 1; first = None
        # sum_{N>Nmax} N^2 PQ z^j sqrt(pi/(2z)) e^{-z + 1/(2z)},  z = c N n; terms decrease geometrically once z > j + 2
        # bound: term_{N+1}/term_N = ((N+1)/N)^{j+2} e^{-cn} * sqrt(N/(N+1)) <= rho
        zN = c*N*n
        t = arb(N)**2*PQ*zN**j*(arb.pi()/(2*zN)).sqrt()*(-zN + 1/(2*zN)).exp()
        rho = (arb(N + 1)/N)**(j + 2)*(-c*n).exp()
        assert rho.upper() < 1, ('tail ratio not < 1', n, j, Nmax)
        tail = t/(1 - rho)
        val = pref*sums[j]
        tb = (pref*tail).upper()
        out.append(val + arb(0, tb))
    return out

def g_basis(n, D, **kw):
    """g_k(log n) = (-1)^k h_{2k}(n), k = 0..D (the function with transform Xi^2 t^{2k})."""
    hv = h_values(n, [2*k for k in range(D + 1)], **kw)
    return [(-1)**k*hv[k] for k in range(D + 1)]
