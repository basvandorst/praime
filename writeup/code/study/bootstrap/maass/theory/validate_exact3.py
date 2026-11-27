# validate_exact3.py: recomputes and asserts every number of exact3.md.  Self-contained (python-flint arb, mpmath, sympy);
# reads only the polynomials B stored in cert_d<deg>.json next to it.
# Run: (ulimit -v 1200000; timeout 900 nice .venv/bin/python -I validate_exact3.py [quick])
#
# Objects.  chi = chi_{-3};  Lambda(s) = 3^{s/2} Gamma_R(s+1) L(s,chi);  Xi(t) = Lambda(1/2+it);  g^ = Xi^2 B, B(t) = sum_k b_k t^{2k}.
#   g(w) = sum_k b_k (-1)^k (phi*phi)^{(2k)}(w),  (phi*phi)(w) = (4/3) e^{3w/2} sum_N chi(N) d(N) N K0(c N e^w),  c = 2 pi/3.
#   P_p(x) = sum_m x^m log p p^{-m/2} g(m log p);  Defect = sum_p [P_p(chi(p)) - min_{[-1,1]} P_p];  bound q >= 3 exp(-2 Defect/g(0)).
import sys; sys.dont_write_bytecode = True
import json, math, time, os
from fractions import Fraction
from decimal import Decimal
from flint import arb, ctx
import mpmath as mp
import sympy as sp

HERE = '/home/bas/praime/writeup/code/study/bootstrap/maass/theory'
QUICK = len(sys.argv) > 1 and sys.argv[1] == 'quick'
T0 = time.time()
NCHK = [0, 0]
def check(name, cond, info=''):
    NCHK[0] += 1
    if not cond: NCHK[1] += 1
    print('%s  %s  %s' % ('ok  ' if cond else 'FAIL', name, info), flush=True)

# ---------------------------------------------------------------- arithmetic
def chi(n):
    r = n % 3
    return 0 if r == 0 else (1 if r == 1 else -1)
def ndiv(n):
    c, d = 0, 1
    while d*d <= n:
        if n % d == 0: c += 1 if d*d == n else 2
        d += 1
    return c
def isprime(n): return n > 1 and all(n % d for d in range(2, int(n**0.5) + 1))
def ppow(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = 0
            while n % p == 0: n //= p; m += 1
            return (p, m) if n == 1 else None
    return None

# ---------------------------------------------------------------- K0, K1 in ball arithmetic
_KC = {}
def _asym(z, nu, prec):
    """K_nu(z), z > 0, from the asymptotic series (DLMF 10.40.2, a_k from 10.17.1) with the bound of DLMF 10.40(ii)
    (read 27 Nov 2026): truncated after l >= max(|nu| - 1/2, 1) terms, the remainder is at most the first neglected term
    in absolute value.  Returns None when the terms stop decreasing before reaching 2^(-prec-8)."""
    mu = 4*nu*nu; term = arb(1); s = arb(0); k = 0; target = arb(2)**(-prec - 8)
    while True:
        s += term; k += 1
        nt = term*(mu - (2*k - 1)**2)/(8*k*z)
        if abs(nt).upper() < target:
            return (arb.pi()/(2*z)).sqrt()*(-z).exp()*(s + arb(0, abs(nt).upper()))
        if abs(nt).upper() > abs(term).upper() and k > 1: return None
        term = nt
def K01(M):
    """K0, K1 at z = 2 pi M/3 (M integer > 0), memoised.  Large z: the asymptotic series with its rigorous remainder.
    Otherwise arb's bessel_k (rigorous balls), with z and K recomputed at a raised precision (arb's series loses about
    2 z log2(e) bits) until the relative radius is below 2^(4-prec)."""
    prec = ctx.prec
    if (M, prec) in _KC: return _KC[(M, prec)]
    z = 2*arb.pi()*M/3
    a0 = _asym(z, 0, prec); a1 = _asym(z, 1, prec)
    if a0 is not None and a1 is not None:
        _KC[(M, prec)] = (a0, a1); return a0, a1
    extra = int(3*2.0944*M) + 60
    while True:
        old = ctx.prec; ctx.prec = prec + extra
        try:
            zz = 2*arb.pi()*M/3; k0 = zz.bessel_k(0); k1 = zz.bessel_k(1)
        finally:
            ctx.prec = old
        if all(k.rad() <= abs(k.mid())*arb(2)**(-prec + 4) for k in (k0, k1)):
            _KC[(M, prec)] = (+k0, +k1); return _KC[(M, prec)]
        extra *= 2; assert extra < 40000

# ---------------------------------------------------------------- the closed form and its derivatives
_PQ = {}
def pq(J):
    """p_j, q_j: D^j[z^{3/2} K0] = z^{3/2}(p_j K0 + q_j K1), D = z d/dz.  Exact rational coefficient lists."""
    if J in _PQ: return _PQ[J]
    p, q = [[Fraction(1)]], [[Fraction(0)]]
    for j in range(J):
        a, b = p[-1], q[-1]; L = max(len(a), len(b)) + 1
        a = a + [Fraction(0)]*(L - len(a)); b = b + [Fraction(0)]*(L - len(b))
        na = [Fraction(0)]*L; nb = [Fraction(0)]*L
        for i in range(L):
            na[i] += Fraction(3, 2)*a[i] + i*a[i]; nb[i] += Fraction(1, 2)*b[i] + i*b[i]
            if i + 1 < L: na[i + 1] -= b[i]; nb[i + 1] -= a[i]
        while len(na) > 1 and na[-1] == 0: na.pop()
        while len(nb) > 1 and nb[-1] == 0: nb.pop()
        p.append(na); q.append(nb)
    _PQ[J] = (p, q); return p, q
def farb(fr): return arb(fr.numerator)/fr.denominator

def g_of_n(n, bcoef):
    """g(log n) for g = sum_k b_k (-1)^k (phi*phi)^{(2k)}, b_k exact rationals; ball, with the tail N > Nmax in the radius.
    Tail: |p_j(z)|, |q_j(z)| <= (sum of |coefficients|) z^j for z >= 1; K_nu(z) <= sqrt(pi/(2z)) e^{-z + 1/(2z)}; d(N) <= N."""
    D = len(bcoef) - 1; J = 2*D; p, q = pq(J); c = 2*arb.pi()/3
    # combined polynomials P(z) = sum_k b_k (-1)^k p_{2k}(z), Q likewise
    P = [Fraction(0)]*(J + 2); Q = [Fraction(0)]*(J + 2); PQbar = Fraction(0)
    for k, bk in enumerate(bcoef):
        s = (-1)**k*bk
        for i, v in enumerate(p[2*k]): P[i] += s*v
        for i, v in enumerate(q[2*k]): Q[i] += s*v
        PQbar += abs(bk)*(sum(abs(v) for v in p[2*k]) + sum(abs(v) for v in q[2*k]))
    Pa = [farb(v) for v in P]; Qa = [farb(v) for v in Q]
    Nmax = 1
    while True:
        z = 2*math.pi/3*Nmax*n
        if z > 2*J + 10 and -z + J*math.log(z) + 3*math.log(Nmax + 1) + math.log(float(PQbar) + 1) < -(ctx.prec + 20)*math.log(2): break
        Nmax += 1
    tot = arb(0)
    for N in range(1, Nmax + 1):
        ch = chi(N)
        if ch == 0: continue
        z = c*N*n; k0, k1 = K01(N*n)
        pv = arb(0); qv = arb(0)
        for i in range(J + 1, -1, -1): pv = pv*z + Pa[i]; qv = qv*z + Qa[i]
        tot += ch*ndiv(N)*N*(pv*k0 + qv*k1)
    N = Nmax + 1; zN = c*N*n
    t = arb(N)**2*farb(PQbar)*zN**J*(arb.pi()/(2*zN)).sqrt()*(-zN + 1/(2*zN)).exp()
    rho = (arb(N + 1)/N)**(J + 2)*(-c*n).exp(); assert rho.upper() < 1
    pref = arb(4)/3*arb(n)*arb(n).sqrt()
    return pref*tot + arb(0, (pref*t/(1 - rho)).upper())

def envelope_tail(bcoef, NC):
    """Rigorous upper bound for 2 sum_{n > NC} log(n) n^{-1/2} |g(log n)| (covers every prime power above NC)."""
    D = len(bcoef) - 1; J = 2*D; p, q = pq(J); c = 2*arb.pi()/3
    PQbar = sum(abs(bk)*(sum(abs(v) for v in p[2*k]) + sum(abs(v) for v in q[2*k])) for k, bk in enumerate(bcoef))
    def E(y):   # bound for |g(log y)|, y >= NC: N-sum dominated geometrically by N = 1
        z = c*y
        t1 = farb(PQbar)*z**J*(arb.pi()/(2*z)).sqrt()*(-z + 1/(2*z)).exp()
        rho = arb(2)**(J + 2)*(-c*y).exp(); assert rho.upper() < 1
        return arb(4)/3*y*arb(y).sqrt()*t1/(1 - rho)
    y = NC + 1
    term = arb(y).log()/arb(y).sqrt()*E(y)
    # ratio of successive terms for n >= y: ((n+1)/n)^{J+3/2} e^{-c} (1 + 1/(n log n)) <= rho2
    rho2 = (arb(y + 1)/y)**(J + 2)*(-c).exp()*(1 + 1/(arb(y)*arb(y).log())); assert rho2.upper() < 1
    return (2*term/(1 - rho2)).upper()

# ---------------------------------------------------------------- rigorous minimum of a polynomial on [-1, 1]
def peval(a, x):   # a[0] + a[1] x + ... (arb coefficients), x an arb (point or interval)
    v = arb(0)
    for ai in reversed(a): v = v*x + ai
    return v
def pder(a): return [i*a[i] for i in range(1, len(a))]
def poly_min_lower(a, maxint=200000):
    """Lower bound of min_{[-1,1]} sum a_i x^i (a_0 = 0 allowed), by branch and bound with derivative signs and a
    convex quadratic Taylor model (P(x) >= P(c) + P'(c)(x-c) + m (x-c)^2/2, m = lower bound of P'' on the interval)."""
    a1 = pder(a); a2 = pder(a1)
    stack = [(arb(-1), arb(1))]; best_upper = min(peval(a, arb(-1)).upper(), peval(a, arb(1)).upper())
    lower = None; count = 0; tol = arb(2)**-200; wmin = arb(2)**-60
    while stack:
        l, r = stack.pop(); count += 1; assert count < maxint
        I = arb((l + r)/2, ((r - l)/2).upper())
        d = peval(a1, I)
        if d > 0:      # increasing on [l, r]: the minimum is P(l)
            v = peval(a, l); lb = v.lower(); best_upper = min(best_upper, v.upper())
        elif d < 0:    # decreasing: the minimum is P(r)
            v = peval(a, r); lb = v.lower(); best_upper = min(best_upper, v.upper())
        else:          # Taylor at the centre: P(x) >= P(c) + P'(c)(x-c) + m (x-c)^2/2 on [l, r], m <= P'' there
            c = (l + r)/2; rad = (r - l)/2
            Pc = peval(a, c); dc = peval(a1, c); m = peval(a2, I).lower()
            best_upper = min(best_upper, Pc.upper())
            if m > 0 and abs(dc/m).lower() >= rad:
                lb = min((Pc + s*dc*rad + m*rad*rad/2).lower() for s in (-1, 1))   # model minimiser surely outside: an end
            elif m > 0:
                lb = (Pc - dc*dc/(2*m)).lower()                 # global minimum of the convex model (always valid)
            else:
                lb = (Pc - abs(dc)*rad + m*rad*rad/2).lower()
            if lb < best_upper - tol and (r - l) > wmin:
                stack.append((l, c)); stack.append((c, r)); continue
        lower = lb if lower is None else min(lower, lb)
    return lower, best_upper, count

# ---------------------------------------------------------------- the certificate for one polynomial B
def certify(bcoef, NC=400, prec=600, verbose=True):
    old = ctx.prec; ctx.prec = prec
    try:
        res = {}
        # 1. B > 0 on R: q(s) = sum b_k s^k has no root on [0, oo) and q(0) > 0 (exact Sturm count, sympy)
        s = sp.symbols('s'); qpoly = sp.Poly([sp.Rational(b.numerator, b.denominator) for b in reversed(bcoef)], s)
        nroots = qpoly.count_roots(0, None) if qpoly.degree() > 0 else 0
        res['Bpos'] = (nroots == 0 and bcoef[0] > 0)
        # 2. g(0) and g at prime powers n <= NC
        g0 = g_of_n(1, bcoef); res['g0'] = g0
        gv = {}
        for n in range(2, NC + 1):
            if ppow(n): gv[n] = g_of_n(n, bcoef)
        # 3. per-prime defects (p <= NC; powers above NC go to the tail)
        tot = arb(0); per = {}
        for p in [p for p in range(2, NC + 1) if isprime(p)]:
            a = [arb(0)]; m = 1
            while p**m <= NC:
                a.append(arb(p).log()*arb(p)**(-arb(m)/2)*gv[p**m]); m += 1
            x = chi(p)
            Pchi = peval(a, arb(x))
            lo, up, cnt = poly_min_lower(a)
            dp = Pchi.upper() - lo        # upper bound of P(chi) - min
            per[p] = (arb(dp), Pchi - arb(up))
            tot += arb(dp)
        tail = envelope_tail(bcoef, NC)
        Dup = (tot + arb(tail)).upper()
        res['Defect_upper'] = Dup; res['tail'] = tail; res['per'] = per; res['gv'] = gv
        ratio_up = (arb(Dup)/g0).upper()
        res['ratio_upper'] = ratio_up
        res['bound_lower'] = (3*(-2*arb(ratio_up)).exp()).lower()
        # a lower estimate of the defect (to show the bound is not loose): sum over p of max(0, P(chi) - upper(min))
        return res
    finally:
        ctx.prec = old

def load_b(fn):
    d = json.load(open(fn))
    return [Fraction(x) for x in d['b_exact']], d

# ---------------------------------------------------------------- A. checks of the closed form (independent: mpmath L-function)
def sec_A():
    print('\n== A. closed form of phi*phi and the identity (independent routes)', flush=True)
    mp.mp.dps = 30; ctx.prec = 200
    def Lam(s): return mp.power(3, s/2)*mp.power(mp.pi, -(s + 1)/2)*mp.gamma((s + 1)/2)*mp.dirichlet(s, [0, 1, -1])
    def theta1(x, K=60): return mp.fsum(n*chi(n)*mp.exp(-mp.pi*n*n*x/3) for n in range(1, K))
    # A1. theta transformation theta1(1/x) = x^{3/2} theta1(x) (root number +1), at x = 0.6, 1.7
    for x in [mp.mpf('0.6'), mp.mpf('1.7')]:
        lhs = theta1(1/x, 200); rhs = x**1.5*theta1(x, 200)
        check('theta transformation at x=%s' % x, abs(lhs - rhs) < mp.mpf(10)**-25, '%s vs %s' % (mp.nstr(lhs, 18), mp.nstr(rhs, 18)))
    # A2. Xi(t) = int phi(u) cos(tu) du (phi from theta1) against Lambda(1/2+it) from mpmath's Dirichlet L, t = 0 and 6.3
    phi = lambda u: 2/mp.sqrt(3)*mp.exp(1.5*u)*theta1(mp.exp(2*u)) if u >= 0 else 2/mp.sqrt(3)*mp.exp(-1.5*u)*theta1(mp.exp(-2*u))
    for t in [mp.mpf(0), mp.mpf('6.3')]:
        I = 2*mp.quad(lambda u: phi(u)*mp.cos(t*u), [0, 0.5, 1, 1.5, 2.5])
        L = Lam(mp.mpf('0.5') + 1j*t)
        check('Mellin: Xi(%s) = int phi cos' % t, abs(I - L) < mp.mpf(10)**-20, mp.nstr(I, 18))
    # A3. (phi*phi)(log 2) by quadrature against the Bessel closed form
    for n in [1, 2, 3]:
        I = mp.quad(lambda v: phi(v + mp.log(n))*phi(v), [-3, -1.5, -0.5, 0, 0.5, 1.5, 3])
        cf = g_of_n(n, [Fraction(1)])
        check('closed form (phi*phi)(log %d) vs quadrature, 20 digits' % n, abs(I - mp.mpf(cf.mid().str(30, radius=False))) < mp.mpf(10)**-20*abs(I), mp.nstr(I, 20))
    old = mp.mp.dps; mp.mp.dps = 130; ctx.prec = 450
    worst = 0
    for M in [1, 7, 60, 100, 150, 400]:
        k0, k1 = K01(M); z = 2*mp.pi*M/3
        for k, nu in ((k0, 0), (k1, 1)):
            ref = mp.besselk(nu, z); worst = max(worst, float(abs(mp.mpf(k.mid().str(125, radius=False)) - ref)/ref))
    mp.mp.dps = old; ctx.prec = 200
    check('K0, K1 (own: arb / asymptotic with DLMF bound) vs mpmath besselk at 130 digits, z = 2.09 .. 838', worst < 1e-115, '%.1e' % worst)
    # A4. the identity 0 = (log 3) g(0) + A_1(g) - 2 sum chi(n) Lambda(n) n^{-1/2} g(log n), with g(0) and A_1 from the
    #     Xi side (mpmath L-function, Gauss-Legendre) and the prime sum from the closed form.  B = 1 and B = 1 + t^2/7.
    cache = {}
    def Xi2(t):
        if t not in cache: cache[t] = Lam(mp.mpf('0.5') + 1j*t).real**2
        return cache[t]
    nodes = [0, 2, 5, 10, 16, 24, 34, 46, 60, 80]
    for k in range(3):
        I = mp.quad(lambda t: Xi2(t)*t**(2*k), nodes)/mp.pi
        e = [Fraction(0)]*k + [Fraction(1)]
        c = mp.mpf(g_of_n(1, e).mid().str(30, radius=False))
        check('g_%d(0) = (1/2pi) int Xi^2 t^%d dt (L-function side) vs closed form' % (k, 2*k), abs(I - c) < mp.mpf(10)**-20*abs(c), mp.nstr(c, 20))
    for name, bc in [('B=1', [Fraction(1)]), ('B=1+t^2/7', [Fraction(1), Fraction(1, 7)])]:
        Bf = lambda t: sum(mp.mpf(b.numerator)/b.denominator*t**(2*k) for k, b in enumerate(bc))
        F = lambda t: Xi2(t)*Bf(t)
        Om = lambda t: -mp.log(mp.pi) + mp.re(mp.digamma(mp.mpf('0.75') + 0.5j*t))   # 2 Re psi_R(3/2 + it)
        nodes = [0, 2, 5, 10, 16, 24, 34, 46, 60, 80]
        g0 = mp.quad(F, nodes)/mp.pi
        A1 = mp.quad(lambda t: F(t)*Om(t), nodes)/mp.pi
        S = mp.mpf(0)
        for n in range(2, 120):
            pm = ppow(n)
            if pm and chi(n) != 0:
                S += chi(n)*mp.log(pm[0])/mp.sqrt(n)*mp.mpf(g_of_n(n, bc).mid().str(30, radius=False))
        g0c = mp.mpf(g_of_n(1, bc).mid().str(30, radius=False))
        idv = mp.log(3)*g0 + A1 - 2*S
        check('g(0) closed form vs Xi side (%s)' % name, abs(g0 - g0c) < mp.mpf(10)**-20, mp.nstr(g0c, 20))
        check('identity (log 3) g(0) + A_1(g) - 2 sum = 0 (%s)' % name, abs(idv) < mp.mpf(10)**-18 * g0, 'value %s, g(0) = %s' % (mp.nstr(idv, 3), mp.nstr(g0, 12)))

# ---------------------------------------------------------------- B. certificates per degree
# The table of exact3.md, section 5: degree -> (Defect/g(0) upper bound, lower estimate, 3 - (certified lower bound for q*),
# the prime carrying the largest share, signs of the coefficients b_k of B).  Upper numbers are rounded up, lower down.
CLAIMS = {
    0:  ('4.532e-2', '4.531e-2', '0.2600', 3, '+'),
    2:  ('6.131e-3', '6.130e-3', '0.03656', 5, '++'),
    4:  ('4.471e-4', '4.470e-4', '2.681e-3', 7, '+++'),
    6:  ('1.381e-5', '1.380e-5', '8.284e-5', 3, '++++'),
    8:  ('5.137e-7', '5.136e-7', '3.082e-6', 13, '+++++'),
    10: ('2.752e-9', '2.751e-9', '1.651e-8', 17, '++++++'),
    12: ('7.480e-11', '7.479e-11', '4.488e-10', 19, '+++++++'),
    14: ('5.253e-13', '5.252e-13', '3.152e-12', 23, '++++++++'),
    16: ('1.511e-17', '1.510e-17', '9.062e-17', 31, '+++++++++'),
    18: ('1.510e-17', '1.509e-17', '9.060e-17', 31, '++++++++++'),
    20: ('1.510e-17', '1.509e-17', '9.060e-17', 31, '+++++++++++'),
    22: ('6.350e-23', '6.349e-23', '3.810e-22', 41, '+++++++++-++'),
}
def sec_B(degs):
    print('\n== B. certificates', flush=True)
    rows = {}
    for deg in degs:
        bc, meta = load_b(os.path.join(HERE, 'cert_d%d.json' % deg))
        t = time.time()
        r = certify(bc, NC=meta.get('NC', 400), prec=meta.get('prec', 600))
        g0 = r['g0']; ratio = arb(r['ratio_upper'])
        low = (sum((max(arb(0), w) for (v, w) in r['per'].values()), arb(0)) - arb(r['tail']))/g0   # lower bound: minus the tail
        gap = 3 - arb(r['bound_lower'])                                            # upper bound for 3 - q*
        top = sorted(r['per'].items(), key=lambda kv: -float(kv[1][0].mid()))
        signs = ''.join('+' if x > 0 else '-' for x in bc)
        print('deg %2d: B>0 %s, g(0) = %s, Defect/g(0) in [%s, %s], 3 - q* <= %s, tail %s [%.0f s]; carriers %s' % (
            deg, r['Bpos'], g0.str(12, radius=False), low.lower().str(6, radius=False), ratio.str(6, radius=False),
            gap.upper().str(6, radius=False), arb(r['tail']).str(3, radius=False), time.time() - t,
            ', '.join('%d: %s' % (p, (v/g0).str(3, radius=False)) for p, (v, w) in top[:3])), flush=True)
        check('deg %d: B > 0 on R (Sturm count 0 on [0, oo) in s = t^2)' % deg, r['Bpos'])
        if deg in CLAIMS:
            up, lo, gp, car, sg = CLAIMS[deg]
            check('deg %d: Defect/g(0) <= %s' % (deg, up), ratio < arb(up))
            check('deg %d: Defect/g(0) >= %s (the bound is not loose)' % (deg, lo), low > arb(lo))
            check('deg %d: q* >= 3 - %s' % (deg, gp), gap < arb(gp))
            check('deg %d: largest share at p = %d' % (deg, car), top[0][0] == car)
            check('deg %d: signs of b_k %s' % (deg, sg), signs == sg)
        rows[deg] = r
    return rows

# ---------------------------------------------------------------- C. structural statements of exact3.md, section 6
DECAY = {0: -1.34, 2: -2.21, 4: -3.35, 6: -4.86, 8: -6.29, 10: -8.56, 12: -10.13, 14: -12.28, 16: -16.82, 18: -16.82, 20: -16.82, 22: -22.20}
def sec_C(rows):
    print('\n== C. structure', flush=True)
    for deg, r in rows.items():
        if deg in DECAY:
            v = float(mp.log10(mp.mpf(arb(r['ratio_upper']).mid().str(20, radius=False))))
            check('deg %d: log10(Defect/g(0)) = %.2f' % (deg, DECAY[deg]), abs(v - DECAY[deg]) < 0.006, '%.4f' % v)
        check('deg %d: tail above 400 < 1e-297' % deg, arb(r['tail']) < arb('1e-297'))
        check('deg %d: Defect_2/g(0) < 1e-100 (p = 2 carries no cost)' % deg, (r['per'][2][0]/r['g0']) < arb('1e-100'))
        if deg >= 2:
            check('deg %d: g(log 2) > 0 > g(log 4)' % deg, r['gv'][2] > 0 and r['gv'][4] < 0)
    # eventual sign of g is (-1)^D sign(b_D) = (-1)^D; the largest share outside p = 3 falls on the class chi(p) = (-1)^D,
    # except at degrees 18 and 20, where the top coefficients are ~0 and the degree-16 pattern persists
    for deg, r in rows.items():
        D = deg//2
        p = max((q for q in r['per'] if q != 3), key=lambda q: float(r['per'][q][0].mid()))
        if deg in (18, 20):
            check('deg %d: largest share outside 3 at p = %d with chi(p) = +1 = (-1)^8 (degree-16 pattern)' % (deg, p), chi(p) == 1)
        else:
            check('deg %d: largest share outside 3 at p = %d with chi(p) = (-1)^D = %+d' % (deg, p, (-1)**D), chi(p) == (-1)**D)
    x0 = float((rows[0]['gv'][4]/rows[0]['g0']).mid()); x2 = float((rows[0]['gv'][2]/rows[0]['g0']).mid())/2**0.5
    check('deg 0: g(log 4)/g(0) = +0.0118 and g(log 2)/(sqrt2 g(0)) = 0.26', abs(x0 - 0.0118) < 0.00005 and abs(x2 - 0.26) < 0.005, '%.5f %.4f' % (x0, x2))
    r = rows[16]; ctx.prec = 600
    want = {2: (0.219, 0.0006), 3: (9e-16, 0.6e-16), 4: (-0.0036, 0.00005), 5: (9e-15, 0.6e-15), 7: (-6e-24, 0.6e-24)}
    for n, (v, tol) in want.items():
        x = float((r['gv'][n]/r['g0']).mid())
        check('deg 16: g(log %d)/g(0) = %g' % (n, v), abs(x - v) < tol, '%.4g' % x)
    for deg, k, lim in [(18, 9, '1e-24'), (20, 9, '1e-24'), (20, 10, '1e-40')]:
        bc, meta = load_b(os.path.join(HERE, 'cert_d%d.json' % deg))
        check('deg %d: b_%d/b_0 in (0, %s)' % (deg, k, lim), 0 < bc[k]/bc[0] < Fraction(Decimal(lim)))

# ---------------------------------------------------------------- D. the repo's LP dual certificate at support 24 (read-only data)
def sec_D():
    print('\n== D. the LP dual of maass/dual.py at support 24: the local polynomial at p = 2', flush=True)
    fn = '/home/bas/praime/writeup/data/study/maass/dual_odd_out.txt'
    G = {}
    for line in open(fn):
        w = line.split()
        if len(w) >= 3 and w[0].isdigit():
            G[int(w[0])] = float(w[2])
    # P_2(x) up to the factor log 2 / G(1): sum_m x^m 2^{-m/2} G(2^m), m = 1..4 (2^m <= 24)
    c = [2**(-m/2)*G[2**m] for m in range(1, 5)]
    slope = sum(m*(-1)**(m - 1)*c[m - 1] for m in range(1, 5))         # P_2'(-1)
    curv = sum(m*(m - 1)*(-1)**(m - 2)*c[m - 1] for m in range(2, 5))  # P_2''(-1)
    print('   G(2), G(4), G(8), G(16) = %s; P_2\'(-1) = %.2e, P_2\'\'(-1) = %.2e (units of G(1))' % ([G[2**m] for m in range(1, 5)], slope, curv))
    check('dual at p=2: P_2\'(-1) = -4.0e-5 (minimum just inside x = -1, defect second order)', abs(slope + 4.03e-5) < 0.2e-5 and curv > 0)
    # the per-power sign test of dual.py flags n = 4 and 16; with G(log 4) > 0 the minimum of P_2 can still sit at x = -1:
    check('dual at p=2: G(4) > 0 and G(16) > 0 (flagged by dual.py) while G(2)/sqrt2 - G(4) is within 1.2e-3 of 0',
          G[4] > 0 and G[16] > 0 and abs(G[2]/2**0.5 - G[4]) < 1.2e-3)

if __name__ == '__main__':
    if not QUICK: sec_A()
    degs = sorted(int(f[6:-5]) for f in os.listdir(HERE) if f.startswith('cert_d') and f.endswith('.json'))
    rows = sec_B(degs)
    sec_C(rows)
    sec_D()
    print('\n%d checks, %d failed   [%.0f s]' % (NCHK[0], NCHK[1], time.time() - T0))
