# validate_exact12.py: recomputes and asserts every number of exact12.md (weight problem at level 1, one form L(f)).
# Self-contained (python-flint arb/acb, mpmath, sympy, cypari2 for tau(n) and an independent Lambda(s, Delta));
# reads only the polynomials B stored in cert12_d<deg>.json.
# Run: (ulimit -v 1200000; timeout 900 nice .venv/bin/python -I validate_exact12.py [quick])
#
# Lambda(s) = Gamma_C(s + 11/2) L(s, Delta) (analytic normalisation), Xi(t) = Lambda(1/2+it), g^ = Xi^2 B.
#   g(w) = sum_k b_k (-1)^k (phi*phi)^{(2k)}(w),  (phi*phi)(w) = 8 e^{6w} sum_N (tau*tau)(N) K0(4 pi sqrt(N) e^{w/2}).
#   A^(k)(g) = (1/pi) int g^(t) [-log 2 pi + Re psi(k/2 + it)] dt  (= A_{(k-1)/2} + A_{(k+1)/2} of maass.md section 2).
#   P_p(x) = sum_m 2 T_m(x) log p p^{-m/2} g(m log p), x = cos theta;  x_p = tau(p)/(2 p^{11/2});
#   Defect = sum_p [P_p(x_p) - min_{[-1,1]} P_p];  A'(12) = int_0^oo u e^{-6u} g(u)/(1 - e^{-u}) du;
#   no feasible weight below 12 - eps, eps = 2 Defect/(A'(12) - delta M2)  (exact12.md, Theorem 4).
import sys; sys.dont_write_bytecode = True
import json, math, time, os
from fractions import Fraction
from flint import arb, acb, ctx
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

# ---------------------------------------------------------------- tau, tau*tau
NT = 8000
import cypari2
PARI = cypari2.Pari()
TAU = [0] + [int(PARI.ramanujantau(n)) for n in range(1, NT + 1)]
TT = [0]*(NT + 1)
for a in range(1, NT + 1):
    for b in range(1, NT//a + 1): TT[a*b] += TAU[a]*TAU[b]
def isprime(n): return n > 1 and all(n % d for d in range(2, int(n**0.5) + 1))
def ppow(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = 0
            while n % p == 0: n //= p; m += 1
            return (p, m) if n == 1 else None
    return None

# ---------------------------------------------------------------- K0, K1 at z = 4 pi sqrt(M)
def _asym(z, nu, prec):
    """DLMF 10.40.2 with the bound of 10.40(ii) (x > 0, l >= max(|nu| - 1/2, 1) terms: remainder <= first neglected term)."""
    mu = 4*nu*nu; term = arb(1); s = arb(0); k = 0; target = arb(2)**(-prec - 8)
    while True:
        s += term; k += 1
        nt = term*(mu - (2*k - 1)**2)/(8*k*z)
        if abs(nt).upper() < target: return (arb.pi()/(2*z)).sqrt()*(-z).exp()*(s + arb(0, abs(nt).upper()))
        if abs(nt).upper() > abs(term).upper() and k > 1: return None
        term = nt
_KC = {}
def K01(M):
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

# ---------------------------------------------------------------- closed form of g and its tails
_PQ = {}
def pq(J, a=12):
    """D^j[z^a K0] = z^a (p_j K0 + q_j K1), D = z d/dz: p' = a p + z p' - z q, q' = (a-1) q + z q' - z p."""
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
_COMB = {}
def combined(bcoef):
    """P_B = sum_k b_k (-1)^k 4^{-k} p_{2k}, Q_B likewise (the factor 2^{-j} of d/dw = (1/2) z d/dz), and sum of |coefficients|."""
    key = tuple(bcoef)
    if key in _COMB: return _COMB[key]
    D = len(bcoef) - 1; J = 2*D; p, q = pq(J)
    P = [Fraction(0)]*(J + 2); Q = [Fraction(0)]*(J + 2)
    for k, bk in enumerate(bcoef):
        s = (-1)**k*bk/4**k
        for i, v in enumerate(p[2*k]): P[i] += s*v
        for i, v in enumerate(q[2*k]): Q[i] += s*v
    bar = sum(abs(v) for v in P) + sum(abs(v) for v in Q)
    _COMB[key] = ([farb(v) for v in P], [farb(v) for v in Q], farb(bar), J)
    return _COMB[key]
def gamma_upper_bound(s, X):
    assert X > s - 1
    return X**(s - 1)*(-X).exp()*X/(X - s + 1)
def powsum_bound(A, b, c, x0):
    """Upper bound for sum_{N > x0} A N^b e^{-c sqrt N} (A, b, c > 0), valid when c sqrt(x0) > 2b (f decreasing there):
    <= int_{x0}^oo A x^b e^{-c sqrt x} dx = 2 A c^{-(2b+2)} Gamma(2b+2, c sqrt x0)."""
    X = c*arb(x0).sqrt(); assert X > 2*b + 1
    return 2*A*c**(-(2*b + 2))*gamma_upper_bound(2*b + 2, X)

def g_of_n(n, bcoef):
    """g(log n), ball.  Tail N > Nmax: |(tau*tau)(N)| <= 8 N^7 (Deligne), |P_B|,|Q_B| <= bar z^J (z >= 1),
    K_nu(z) <= sqrt(pi/(2z)) e^{-z + 1/(2z)}."""
    Pa, Qa, bar, J = combined(bcoef)
    Nmax = 1
    while True:
        z = 4*math.pi*math.sqrt(Nmax*n)
        if z > 2*(7 + J/2) + 20 and -z + J*math.log(z) + 8*math.log(Nmax + 1) + math.log(float(bar.mid()) + 1) < -(ctx.prec + 20)*math.log(2) - 10: break
        Nmax += 1
    assert Nmax <= NT
    tot = arb(0)
    for N in range(1, Nmax + 1):
        if TT[N] == 0: continue
        z = 4*arb.pi()*arb(N*n).sqrt(); k0, k1 = K01(N*n)
        pv = arb(0); qv = arb(0)
        for i in range(J + 1, -1, -1): pv = pv*z + Pa[i]; qv = qv*z + Qa[i]
        tot += TT[N]*(pv*k0 + qv*k1)
    c = 4*arb.pi()*arb(n).sqrt(); b = arb(7) + arb(J)/2 - arb(1)/4
    A = 8*bar*c**J*(arb.pi()/(2*c)).sqrt()*(1/(2*c)).exp()
    t = powsum_bound(A, b, c, Nmax)
    pref = 8*arb(n)**6
    return pref*tot + arb(0, (pref*t).upper())

def envelope_tail(bcoef, NC):
    """Upper bound for 4 sum_{n > NC} log n n^{-1/2} |g(log n)|, the defect of everything above NC (|2 T_m| <= 2 on [-1, 1]).
    With M = nN: log n n^{-1/2} 8 n^6 8 N^7 <= 64 M^7 log M, the number of (n, N) with nN = M is d(M) <= 2 sqrt M, log M <= sqrt M."""
    Pa, Qa, bar, J = combined(bcoef)
    c = 4*arb.pi(); b = arb(8) + arb(J)/2 - arb(1)/4
    A = 128*bar*c**J*(arb.pi()/(2*c)).sqrt()*(1/(2*c*arb(NC).sqrt())).exp()
    return (4*powsum_bound(A, b, c, NC)).upper()

# ---------------------------------------------------------------- rigorous minimum on [-1, 1] (repaired branch m > 0)
def peval(a, x):
    v = arb(0)
    for ai in reversed(a): v = v*x + ai
    return v
def pder(a): return [i*a[i] for i in range(1, len(a))]
def poly_min_lower(a, maxint=200000):
    a1 = pder(a); a2 = pder(a1)
    stack = [(arb(-1), arb(1))]; best_upper = min(peval(a, arb(-1)).upper(), peval(a, arb(1)).upper())
    lower = None; count = 0; tol = arb(2)**-200; wmin = arb(2)**-60
    while stack:
        l, r = stack.pop(); count += 1; assert count < maxint
        I = arb((l + r)/2, ((r - l)/2).upper())
        d = peval(a1, I)
        if d > 0:
            v = peval(a, l); lb = v.lower(); best_upper = min(best_upper, v.upper())
        elif d < 0:
            v = peval(a, r); lb = v.lower(); best_upper = min(best_upper, v.upper())
        else:
            c = (l + r)/2; rad = (r - l)/2
            Pc = peval(a, c); dc = peval(a1, c); m = peval(a2, I).lower()
            best_upper = min(best_upper, Pc.upper())
            if m > 0 and abs(dc/m).lower() >= rad:
                lb = min((Pc + s*dc*rad + m*rad*rad/2).lower() for s in (-1, 1))
            elif m > 0:
                lb = (Pc - dc*dc/(2*m)).lower()
            else:
                lb = (Pc - abs(dc)*rad + m*rad*rad/2).lower()
            if lb < best_upper - tol and (r - l) > wmin:
                stack.append((l, c)); stack.append((c, r)); continue
        lower = lb if lower is None else min(lower, lb)
    return lower, best_upper, count
def cheb(M):
    T = [[1], [0, 1]]
    while len(T) <= M:
        a, b = T[-1], T[-2]; n = [0] + [2*v for v in a]
        for i, v in enumerate(b): n[i] -= v
        T.append(n)
    return T

# ---------------------------------------------------------------- A'(12) = int_0^oo u e^{-6u} g(u)/(1 - e^{-u}) du: a rigorous lower bound
def hurwitz(n, a):
    """zeta(n, a) for integer a >= 1: zeta(n) - sum_{k<a} k^{-n}."""
    return arb(n).zeta() - sum((arb(k)**(-n) for k in range(1, a)), arb(0))
def Aprime12_lower(bcoef, jmax=5):
    """A'(12) = int_0^oo h(s) g(s) ds, h(s) = s e^{-6s}/(1-e^{-s}) >= 0, g(s) = (1/2pi) int g^(t) cos(ts) dt, g^ >= 0.
    cos x >= sum_{j<=2m+1} (-1)^j x^{2j}/(2j)! for all real x, and int_0^oo h(s) s^{2j} ds = (2j+1)! zeta(2j+2, 6), so
    A'(12) >= sum_{j<=2m+1} (-1)^j (2j+1) zeta(2j+2, 6) M_{2j},  M_{2j} = (1/2pi) int g^ t^{2j} dt = g_{B t^{2j}}(0) (closed form).
    Returns the best of these lower bounds over m with 2m+1 <= jmax, and the index used."""
    S = arb(0); best = None; jb = None
    for j in range(jmax + 1):
        Mj = g_of_n(1, [Fraction(0)]*j + list(bcoef))
        S += (-1)**j*(2*j + 1)*hurwitz(2*j + 2, 6)*Mj
        if j % 2 == 1 and (best is None or S.lower() > best):
            best = S.lower(); jb = j
    return best, jb

# ---------------------------------------------------------------- the certificate for one polynomial B
def certify(bcoef, NC=700, prec=600):
    old = ctx.prec; ctx.prec = prec
    try:
        res = {}
        s = sp.symbols('s'); qpoly = sp.Poly([sp.Rational(b.numerator, b.denominator) for b in reversed(bcoef)], s)
        nroots = qpoly.count_roots(0, None) if qpoly.degree() > 0 else 0
        res['Bpos'] = (nroots == 0 and bcoef[0] > 0)
        g0 = g_of_n(1, bcoef); res['g0'] = g0
        gv = {n: g_of_n(n, bcoef) for n in range(2, NC + 1) if ppow(n)}
        res['gv'] = gv
        tot = arb(0); low = arb(0); per = {}
        for p in [p for p in range(2, NC + 1) if isprime(p)]:
            M = 1
            while p**(M + 1) <= NC: M += 1
            T = cheb(M)
            cm = [arb(p).log()*arb(p)**(-arb(m)/2)*gv[p**m] for m in range(1, M + 1)]
            a = [arb(0)] + [sum((2*T[m][k]*cm[m - 1] for m in range(k, M + 1) if k < len(T[m])), arb(0)) for k in range(1, M + 1)]
            xp = arb(TAU[p])/(2*arb(p)**5*arb(p).sqrt())
            Pxp = peval(a, xp)
            lo, up, cnt = poly_min_lower(a)
            dp = arb(Pxp.upper() - lo)
            per[p] = (dp, Pxp - arb(up), a)
            tot += dp; low += max(arb(0), Pxp - arb(up))
        tail = envelope_tail(bcoef, NC)
        Dup = (tot + arb(tail)).upper()
        res.update(Defect_upper=Dup, Defect_lower=(low - arb(tail)).lower(), tail=tail, per=per)
        return res
    finally:
        ctx.prec = old

def Delta_lower(Ms, kappa):
    """Rigorous lower bound for Delta_B(kappa) = A^(12)(g) - A^(kappa)(g) = int_kappa^12 A'(s) ds, 0 < kappa < 12.
    A'(s) = int_0^oo h_s(u) g(u) du with h_s(u) = u e^{-su/2}/(1-e^{-u}) >= 0, and the cosine bound (Aprime12_lower) give
    A'(s) >= sum_{j<=2m+1} (-1)^j (2j+1) zeta(2j+2, s/2) M_{2j}.  Integrated over s in [kappa, 12]
    (d/da zeta(n,a) = -n zeta(n+1,a), zeta(2,a) = psi'(a)):
    Delta_B(kappa) >= 2 M_0 [psi(6) - psi(kappa/2)] + 2 sum_{j=1}^{2m+1} (-1)^j M_{2j} [zeta(2j+1, kappa/2) - zeta(2j+1, 6)].
    Returns the best lower end over the odd truncations 2m+1 < len(Ms)."""
    a = arb(kappa)/2; best = None
    S = 2*Ms[0]*(arb(6).digamma() - a.digamma())
    for j in range(1, len(Ms)):
        S += 2*(-1)**j*Ms[j]*(arb(2*j + 1).zeta(a) - arb(2*j + 1).zeta(arb(6)))
        if j % 2 == 1 and (best is None or S.lower() > best): best = S.lower()
    return arb(best)

def weight_eps(Dup, Ms):
    """Small eps with Delta_lower(12 - eps) > 2 Defect_upper, certified in ball arithmetic: scan eps = 1e-30 * 2^k upward
    to the first success (eps <= 8), then bisect in log eps.  None if no eps <= 8 works."""
    target = 2*arb(Dup)
    ok = lambda e: Delta_lower(Ms, 12 - arb(e)) > target
    lo, hi = None, None; e = 1e-30
    while e <= 8:
        if ok(e): hi = e; break
        lo = e; e *= 2
    if hi is None: return None
    if lo is not None:
        for _ in range(40):
            mid = math.sqrt(lo*hi)
            if ok(mid): hi = mid
            else: lo = mid
    assert ok(hi)
    return arb(hi)

def load_b(fn):
    d = json.load(open(fn)); return [Fraction(x) for x in d['b_exact']], d

# ---------------------------------------------------------------- A. checks of the set-up (independent routes)
def sec_A():
    print('\n== A. Delta: q-expansion, modularity, Mellin, PARI L-function, closed form, Lemma 1, Polya', flush=True)
    mp.mp.dps = 30; ctx.prec = 200
    # A1. tau from PARI against the own q-expansion q prod (1 - q^n)^24, n <= 150
    N = 150; s = [0]*(N + 1); s[0] = 1
    for n in range(1, N + 1):
        for _ in range(24):
            for k in range(N, n - 1, -1): s[k] -= s[k - n]
    check('tau(n), n <= 150: PARI = own q-expansion', [0] + s[:N] == TAU[:N + 1])
    Dl = lambda y: mp.fsum(TAU[n]*mp.exp(-2*mp.pi*n*y) for n in range(1, 400))
    y = mp.mpf('1.3')
    check('Delta(i/y) = y^12 Delta(iy) at y = 1.3', abs(Dl(1/y) - y**12*Dl(y)) < mp.mpf(10)**-25*abs(Dl(1/y)))
    phi = lambda u: 2*mp.exp(6*u)*Dl(mp.exp(u)) if u >= 0 else 2*mp.exp(-6*u)*Dl(mp.exp(-u))
    s8 = mp.mpf(8)
    M = mp.quad(lambda u: phi(u)*mp.exp((s8 - mp.mpf('0.5'))*u), [-3, -1, 0, 1, 3])
    L = mp.fsum(TAU[n]*mp.power(n, -mp.mpf(11)/2 - s8) for n in range(1, NT + 1))
    G = 2*(2*mp.pi)**(-(s8 + mp.mpf(11)/2))*mp.gamma(s8 + mp.mpf(11)/2)*L
    check('Mellin: Lambda(8) = 2 int Delta(iy) y^{8+11/2} dy/y = Gamma_C(8+11/2) L(8)', abs(M - G) < mp.mpf(10)**-20*abs(G), mp.nstr(M, 18))
    # PARI's Lambda (arithmetic normalisation: Lambda_PARI(s) = Gamma_C(s) L_ar(s), so Xi(t) = Lambda_PARI(6 + it))
    PARI.set_real_precision(38)
    Ld = PARI('lfuncreate([n->vector(n,i,ramanujantau(i)),0,[0,1],12,1,1])')   # arithmetic normalisation, Gamma_C(s), s <-> 12 - s
    check('PARI lfuncheckfeq for Delta <= -100 (functional equation to 2^-100)', int(PARI('lfuncheckfeq')(Ld)) <= -100)
    XiP = lambda t: mp.mpf(str(PARI.real(PARI.lfunlambda(Ld, PARI('6 + I*' + mp.nstr(t, 30))))).replace(' E', 'e'))
    for t in [mp.mpf(0), mp.mpf('5.5')]:
        I = 2*mp.quad(lambda u: phi(u)*mp.cos(t*u), [0, 0.5, 1, 2, 3])
        check('Xi(%s): own Mellin = PARI lfunlambda(Delta, 6 + it), 17 digits' % t, abs(I - XiP(t)) < mp.mpf(10)**-17*abs(I), '%s vs %s' % (mp.nstr(I, 20), mp.nstr(XiP(t), 20)))
    for n in [1, 2]:
        I = mp.quad(lambda v: phi(v + mp.log(n))*phi(v), [-3, -1.5, -0.5, 0, 0.5, 1.5, 3])
        cf = g_of_n(n, [Fraction(1)])
        check('closed form (phi*phi)(log %d) vs quadrature, 20 digits' % n, abs(I - mp.mpf(cf.mid().str(30, radius=False))) < mp.mpf(10)**-20*abs(I), mp.nstr(I, 20))
    # Lemma 1: A^(12)(g) on the Xi side (PARI) = 2 sum_p P_p(x_p) (closed form), B = 1 and B = 1 + t^2/50
    cache = {}
    def Xi2(t):
        if t not in cache: cache[t] = XiP(t)**2
        return cache[t]
    nodes = [0, 2, 5, 10, 16, 24, 34, 46, 60]
    for name, bc in [('B=1', [Fraction(1)]), ('B=1+t^2/50', [Fraction(1), Fraction(1, 50)])]:
        Bf = lambda t: sum(mp.mpf(b.numerator)/b.denominator*t**(2*k) for k, b in enumerate(bc))
        A12 = 2*mp.quad(lambda t: Xi2(t)*Bf(t)*(-mp.log(2*mp.pi) + mp.re(mp.digamma(6 + 1j*t))), nodes)/mp.pi
        g0 = 2*mp.quad(lambda t: Xi2(t)*Bf(t), nodes)/(2*mp.pi)
        S = mp.mpf(0)
        for p in [p for p in range(2, 400) if isprime(p)]:
            m = 1; xp = mp.mpf(TAU[p])/(2*mp.power(p, mp.mpf(11)/2))
            while p**m < 400:
                S += mp.chebyt(m, xp)*2*mp.log(p)*mp.power(p, -mp.mpf(m)/2)*mp.mpf(g_of_n(p**m, bc).mid().str(30, radius=False)); m += 1
        g0c = mp.mpf(g_of_n(1, bc).mid().str(30, radius=False))
        check('g(0) closed form vs Xi side, PARI (%s)' % name, abs(g0 - g0c) < mp.mpf(10)**-18*g0, mp.nstr(g0, 15))
        check('Lemma 1: A^(12)(g) = 2 sum_p P_p(x_p) (%s)' % name, abs(A12 - 2*S) < mp.mpf(10)**-16*g0, 'A^(12) = %s, residual %s' % (mp.nstr(A12, 15), mp.nstr(A12 - 2*S, 3)))
    # Polya lemma (proved for x >= 2); numerical sample, and the failure below 1/2 that forces the restriction
    worst = min(mp.re(mp.psi(1, x + 1j*y)) for x in [2, 3, 6] for y in [0.5*k for k in range(0, 400)])
    check('Re psi\'(x+iy) >= 0 at x = 2, 3, 6, y = 0..199.5 (sample of the proved lemma)', worst >= 0, mp.nstr(worst, 3))
    neg = mp.re(mp.psi(1, mp.mpf('0.3') + 50j))
    check('Re psi\'(0.3 + 50i) < 0 (monotonicity fails for small x: the restriction kappa >= 4 is needed)', neg < 0, mp.nstr(neg, 3))

# ---------------------------------------------------------------- B. certificates
CLAIMS = {   # degree: (Defect/g(0) upper, lower, A'(12)/g(0) lower bound, eps upper, largest share, signs of b_k)
    0:  ('0.3584', '0.3583', '0.1600', '4.086', 2, '+'),
    2:  ('0.1828', '0.1827', '0.1526', '2.317', 2, '++'),
    4:  ('0.04740', '0.04739', '0.1324', '0.7236', 5, '+++'),
    6:  ('8.044e-3', '8.043e-3', '0.1244', '0.1293', 7, '++++'),
    8:  ('5.403e-4', '5.402e-4', '0.1230', '8.785e-3', 11, '+++++'),
    10: ('3.848e-5', '3.847e-5', '0.1224', '6.285e-4', 13, '++++++'),
    12: ('2.112e-6', '2.111e-6', '0.1221', '3.459e-5', 19, '+++++++'),
    14: ('1.332e-7', '1.331e-7', '0.1219', '2.185e-6', 19, '++++++++'),
    16: ('2.324e-8', '2.323e-8', '0.1218', '3.813e-7', 19, '+++++++++'),
    18: ('2.249e-9', '2.248e-9', '0.1217', '3.694e-8', 5, '++++++++++'),
}
def sec_B(degs):
    print('\n== B. certificates', flush=True)
    rows = {}
    for deg in degs:
        bc, meta = load_b(os.path.join(HERE, 'cert12_d%d.json' % deg))
        t = time.time()
        r = certify(bc)
        ctx.prec = 600
        Ap, jb = Aprime12_lower(bc); Ap = arb(Ap)
        Ms = [g_of_n(1, [Fraction(0)]*j + list(bc)) for j in range(8)]
        g0 = r['g0']
        eps = weight_eps(r['Defect_upper'], Ms)
        top = sorted(r['per'].items(), key=lambda kv: -float(kv[1][0].mid()))
        signs = ''.join('+' if x > 0 else '-' for x in bc)
        print('deg %2d: B>0 %s, g(0) = %s, Defect/g(0) in [%s, %s], A\'(12)/g(0) >= %s (j <= %d), eps <= %s, tail %s [%.0f s]; largest %s' % (
            deg, r['Bpos'], g0.str(10, radius=False), (arb(r['Defect_lower'])/g0).str(6, radius=False), (arb(r['Defect_upper'])/g0).str(6, radius=False),
            (Ap/g0).lower().str(6, radius=False), jb, arb(eps).str(6, radius=False) if eps is not None else 'none', arb(r['tail']).str(2, radius=False), time.time() - t,
            ', '.join('%d: %s' % (p, (v[0]/g0).str(3, radius=False)) for p, v in top[:3])), flush=True)
        check('deg %d: B > 0 on R (Sturm)' % deg, r['Bpos'])
        if deg in CLAIMS:
            up, lo, ap, ep, car, sg = CLAIMS[deg]
            check('deg %d: Defect/g(0) <= %s' % (deg, up), arb(r['Defect_upper'])/g0 < arb(up))
            check('deg %d: Defect/g(0) >= %s' % (deg, lo), arb(r['Defect_lower'])/g0 > arb(lo))
            check('deg %d: A\'(12)/g(0) >= %s' % (deg, ap), Ap/g0 > arb(ap))
            check('deg %d: eps <= %s' % (deg, ep), eps is not None and arb(eps) < arb(ep))
            check('deg %d: largest share at p = %d' % (deg, car), top[0][0] == car)
            check('deg %d: signs of b_k %s' % (deg, sg), signs == sg)
        rows[deg] = (r, Ap, eps)
    return rows

def sec_C(rows):
    print('\n== C. the theorem', flush=True)
    eps1 = rows[0][2]
    check('B = 1 covers (0, 12 - eps_1) and 12 - eps_1 > 4 (so the monotone range [4, 12] takes over)', 12 - arb(eps1) > 4, (12 - arb(eps1)).str(6, radius=False))
    check('B = 1: g > 0 (phi > 0): g(log n) > 0 at all prime powers n <= 700 (sample of the proved positivity)', all(v > 0 for v in rows[0][0]['gv'].values()))
    # one form pins theta_2: at degree 8 the minimum of P_2 on [-1, 1] lies next to x_2 = tau(2)/(2 2^{11/2})
    a = rows[8][0]['per'][2][2]; xs = [-1 + k/20000 for k in range(40001)]
    af = [float(c.mid()) for c in a]
    vals = [sum(af[i]*x**i for i in range(len(af))) for x in xs]; xm = xs[vals.index(min(vals))]
    x2 = -24/(2*2**5.5)
    check('deg 8: argmin of P_2 = -0.268 (grid 5e-5), x_2 = %.4f' % x2, abs(xm + 0.268) < 0.0006 and abs(x2 + 0.2652) < 0.0001, '%.4f' % xm)
    d2 = rows[8][0]['per'][2][0]/rows[8][0]['g0']
    check('deg 8: Defect_2/g(0) in (1.1e-7, 1.2e-7)', d2 > arb('1.1e-7') and d2 < arb('1.2e-7'), d2.str(4, radius=False))
    best = min((d for d in rows if rows[d][2] is not None), key=lambda d: float(arb(rows[d][2]).mid()))
    print('   best degree %d: no feasible weight below 12 - %s' % (best, arb(rows[best][2]).str(6, radius=False)))

if __name__ == '__main__':
    if not QUICK: sec_A()
    degs = sorted(int(f[8:-5]) for f in os.listdir(HERE) if f.startswith('cert12_d') and f.endswith('.json'))
    rows = sec_B(degs)
    sec_C(rows)
    print('\n%d checks, %d failed   [%.0f s]' % (NCHK[0], NCHK[1], time.time() - T0))
