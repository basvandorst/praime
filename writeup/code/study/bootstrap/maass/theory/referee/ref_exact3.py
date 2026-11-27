# Referee's independent recomputation for exact3.md (own code; mpmath only; no import of the author's scripts).
#
# Route: g(u) = (1/2pi) int ghat(t) e^{-itu} dt with ghat(t) = Xi(t)^2 B(t), Xi(t) = Lambda(1/2+it),
#        Lambda(s) = 3^{s/2} pi^{-(s+1)/2} Gamma((s+1)/2) L(s, chi_{-3}), L from mpmath.dirichlet (Hurwitz zeta).
# Quadrature: trapezoid with step h on the whole line.  By Poisson summation
#        (h/2pi) sum_k ghat(kh) e^{-ikhu} = sum_j g(u + 2 pi j/h),
# so the only errors are aliasing (g(2pi/h - u), u <= log 400 < 6, 2pi/h = 25.1: ~ exp(-2 e^19)) and the cut at |t| = T.
# A_1(g) from its definition in u (Gauss-Legendre panels); prime side sum chi(n) Lambda(n) n^{-1/2} g(log n).
# Min of P_p on [-1, 1]: endpoints + every real critical point (polyroots of P_p', cross-checked by a grid + Newton).
# Negative control: ghat = Xi_{chi_-4}^2 with the chi_{-3} explicit formula; the residual must equal the zero sum
#        2 sum_{gamma > 0} Xi_{-4}(gamma)^2 over the zeros of L(s, chi_{-3}) (found here as sign changes of Xi_{-3}).
# Run: (ulimit -v 1200000; timeout 600 nice /home/bas/praime/.venv/bin/python -I ref_exact3.py)
import sys, json, time
from fractions import Fraction
import mpmath as mp
from mpmath.calculus.quadrature import GaussLegendre

mp.mp.dps = 60
T0 = time.time()
CERT = '/home/bas/praime/writeup/code/study/bootstrap/maass/theory/cert_d%d.json'
NCHK = [0, 0]
def check(name, cond, info=''):
    NCHK[0] += 1
    if not cond: NCHK[1] += 1
    print('%s  %s  %s' % ('ok  ' if cond else 'FAIL', name, info), flush=True)

def chi3(n):
    r = n % 3
    return 0 if r == 0 else (1 if r == 1 else -1)
def is_prime(n): return n > 1 and all(n % d for d in range(2, int(n**0.5) + 1))
def ppow(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = 0
            while n % p == 0: n //= p; m += 1
            return (p, m) if n == 1 else None
    return None

def Lam(s, q=3, chi=(0, 1, -1)):
    return mp.power(q, s/2)*mp.power(mp.pi, -(s + 1)/2)*mp.gamma((s + 1)/2)*mp.dirichlet(s, list(chi))

# ---------------------------------------------------------------- Xi on the grid t = k h
h = mp.mpf(1)/4; K = 960            # T = 240
def xi_grid(q, chi):
    xs, worst = [], mp.mpf(0)
    for k in range(K + 1):
        v = Lam(mp.mpf(1)/2 + 1j*k*h, q, chi)
        xs.append(v.real)
        if abs(v.real) > mp.mpf(10)**-150: worst = max(worst, abs(v.imag)/abs(v.real))
    return xs, worst
XI3, w3 = xi_grid(3, (0, 1, -1))
print('Xi_{-3} on %d nodes, max |Im|/|Re| = %s  [%.0f s]' % (K + 1, mp.nstr(w3, 3), time.time() - T0), flush=True)
check('Xi_{-3}(t) real on the grid (functional equation, root number +1)', w3 < mp.mpf(10)**-45, mp.nstr(w3, 3))

def Bval(b, t): return mp.fsum(mp.mpf(x.numerator)/x.denominator*t**(2*k) for k, x in enumerate(b))
def ghat_grid(b, XI): return [XI[k]**2*Bval(b, k*h) for k in range(K + 1)]
def gfun(gh, u):
    """(h/2pi)[gh_0 + 2 sum_k gh_k cos(k h u)], Chebyshev recurrence for cos(k h u)."""
    c1 = mp.cos(h*u); cm, ck = mp.mpf(1), c1
    s = gh[1]*c1
    for k in range(2, K + 1):
        cm, ck = ck, 2*c1*ck - cm
        s += gh[k]*ck
    return h/(2*mp.pi)*(gh[0] + 2*s)

GL48 = GaussLegendre(mp.mp).calc_nodes(5, mp.mp.prec)
GL96 = GaussLegendre(mp.mp).calc_nodes(6, mp.mp.prec)
PANELS = [0, 0.25, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6.5]
def omega_weights(nodes):
    """A_1(g) = -(gamma + log pi) g(0) + int_0^oo [2 e^{-2u} g(0) - 2 e^{-3u/2} g(u)]/(1 - e^{-2u}) du (g even, real),
    the u-integral on [0, U] by Gauss-Legendre panels, U = 6.5; on u > U the first term exactly (-g(0) log(1 - e^{-2U})),
    the second is below |g(6.5)| ~ exp(-1300).  With g = (h/2pi) sum_k eps_k gh_k cos(khu) (finite sum) this is
    A_1 = (h/2pi) sum_k eps_k gh_k Om[k]; Om[k] is computed once here (same finite sums, reordered)."""
    U = mp.mpf(PANELS[-1])
    Om = [-(mp.euler + mp.log(mp.pi)) - mp.log(1 - mp.exp(-2*U))]*(K + 1)
    for a, b in zip(PANELS[:-1], PANELS[1:]):
        a, b = mp.mpf(a), mp.mpf(b); hw = (b - a)/2; c = (a + b)/2
        for x, w in nodes:
            u = c + hw*x; den = -mp.expm1(-2*u); e2 = mp.exp(-2*u); e3 = mp.exp(-1.5*u)
            c1 = mp.cos(h*u); cm, ck = mp.mpf(1), c1
            W = hw*w*2/den
            Om[0] += W*(e2 - e3)
            for k in range(1, K + 1):
                Om[k] += W*(e2 - e3*ck)
                cm, ck = ck, 2*c1*ck - cm
    return Om
OM48 = omega_weights(GL48); OM96 = omega_weights(GL96)
print('A_1 weights done [%.0f s]' % (time.time() - T0), flush=True)
def A1(gh, Om):
    return h/(2*mp.pi)*(gh[0]*Om[0] + 2*mp.fsum(gh[k]*Om[k] for k in range(1, K + 1)))

def A1_tside(b, chiq=3, chi=(0, 1, -1)):
    """Second route for A_1 (only as a check of the Gauss step): (1/pi) int_0^oo Xi^2 B Omega dt,
    Omega(t) = -log pi + Re psi(3/4 + it/2), at 30 digits with fresh L-values."""
    with mp.workdps(30):
        f = lambda t: Lam(mp.mpf(1)/2 + 1j*t, chiq, chi).real**2*Bval(b, t)*(-mp.log(mp.pi) + mp.re(mp.digamma(mp.mpf(3)/4 + 1j*t/2)))
        I = mp.quad(f, [0, 2, 5, 10, 16, 24, 34, 46, 60, 80, 110, 150])/mp.pi
        g0 = mp.quad(lambda t: Lam(mp.mpf(1)/2 + 1j*t, chiq, chi).real**2*Bval(b, t), [0, 2, 5, 10, 16, 24, 34, 46, 60, 80, 110, 150])/mp.pi
    return I, g0

PP = [n for n in range(2, 401) if ppow(n)]
PRIMES = [p for p in range(2, 401) if is_prime(p)]

def pmin(c):
    """min over [-1,1] of P(x) = sum_{m>=1} c[m-1] x^m: endpoints and real critical points (two methods)."""
    P = lambda x: mp.fsum(cm*x**(m + 1) for m, cm in enumerate(c))
    dP = lambda x: mp.fsum((m + 1)*cm*x**m for m, cm in enumerate(c))
    cand = [mp.mpf(-1), mp.mpf(1)]
    scale = max(abs(x) for x in c)
    cc = [x for x in c]
    while len(cc) > 1 and abs(cc[-1]) < scale*mp.mpf(10)**-50: cc.pop()
    if len(cc) >= 2 and scale > 0:
        der = [(m + 1)*cc[m]/scale for m in range(len(cc))]          # ascending
        try:
            rts = mp.polyroots(list(reversed(der)), maxsteps=400, extraprec=200)
        except mp.libmp.NoConvergence:
            rts = []
        for r in rts:
            if abs(mp.im(r)) < mp.mpf(10)**-25 and -1 < mp.re(r) < 1: cand.append(mp.re(r))
    v1 = min(P(x) for x in cand)
    # grid + Newton (independent of polyroots)
    N = 2000; xs = [mp.mpf(-1) + mp.mpf(2)*i/N for i in range(N + 1)]; vs = [P(x) for x in xs]
    cand2 = [mp.mpf(-1), mp.mpf(1)]
    for i in range(0, N + 1):      # local minima of the samples, endpoints included (an interior minimum can sit
                                   # between -1 and the first grid point: p = 5 at degrees 16 and 22, x* = -0.99992)
        if (i == 0 or vs[i] <= vs[i - 1]) and (i == N or vs[i] <= vs[i + 1]):
            try:
                r = mp.findroot(dP, xs[i])
                if -1 < r < 1: cand2.append(r)
            except Exception:
                cand2.append(xs[i])
    v2 = min(P(x) for x in cand2)
    return v1, v2

def defect(gv, g0):
    per = {}
    for p in PRIMES:
        c = []; m = 1
        while p**m <= 400:
            c.append(mp.log(p)*mp.power(p, -mp.mpf(m)/2)*gv[p**m]); m += 1
        Pchi = mp.fsum(cm*chi3(p)**(m + 1) for m, cm in enumerate(c))
        v1, v2 = pmin(c)
        per[p] = (Pchi - min(v1, v2), abs(v1 - v2))
        if abs(v1 - v2) > mp.mpf(10)**-40*g0:
            print('   p = %d: polyroots min %s, grid+Newton min %s, P(chi) %s, coeffs/g0 %s' % (p, mp.nstr(v1/g0, 12), mp.nstr(v2/g0, 12),
                  mp.nstr(Pchi/g0, 12), [mp.nstr(x/g0, 6) for x in c]), flush=True)
    return per

def run(b, name, XI=XI3, full=True):
    t = time.time()
    gh = ghat_grid(b, XI)
    cut = abs(gh[-1])/max(abs(x) for x in gh)
    g0 = gfun(gh, mp.mpf(0))
    gv = {n: gfun(gh, mp.log(n)) for n in PP}
    S = mp.fsum(chi3(n)*mp.log(ppow(n)[0])/mp.sqrt(n)*gv[n] for n in PP if chi3(n) != 0)
    a48 = A1(gh, OM48); a96 = A1(gh, OM96)
    res = mp.log(3)*g0 + a96 - 2*S
    out = {'g0': g0, 'gv': gv, 'A1': a96, 'A1diff': abs(a48 - a96), 'res': res, 'cut': cut, 'gh': gh}
    if full:
        per = defect(gv, g0); out['per'] = per
        D = mp.fsum(v for v, _ in per.values()); out['D'] = D
        out['ratio'] = D/g0; out['bound'] = 3*mp.exp(-2*D/g0)
        out['top'] = sorted(per, key=lambda p: -per[p][0])[:3]
        out['mindis'] = max(e for _, e in per.values())
    print('%s: g(0) = %s, A1 = %s (48 vs 96 nodes %s), ghat cut %s, identity residual %s [%.0f s]' % (
        name, mp.nstr(g0, 15), mp.nstr(a96, 15), mp.nstr(out['A1diff'], 2), mp.nstr(cut, 2), mp.nstr(res, 3), time.time() - t), flush=True)
    if full:
        print('   Defect/g(0) = %s, 3 - bound = %s, carriers %s, min methods differ by <= %s' % (
            mp.nstr(out['ratio'], 10), mp.nstr(3 - out['bound'], 10),
            ', '.join('%d: %s' % (p, mp.nstr(out['per'][p][0]/g0, 4)) for p in out['top']), mp.nstr(out['mindis'], 2)), flush=True)
    return out

def load(deg):
    d = json.load(open(CERT % deg))
    return [Fraction(x) for x in d['b_exact']]

# ---------------------------------------------------------------- 1. B = 1: closed form (own besselk), quadrature checks
b0 = [Fraction(1)]
r0 = run(b0, 'B=1')
def phiphi(n):   # (phi*phi)(log n) = (4/3) n^{3/2} sum_N chi(N) d(N) N K0(2 pi N n/3), own code
    s = mp.mpf(0)
    for N in range(1, 200):
        if chi3(N) == 0: continue
        dN = sum(1 for d in range(1, N + 1) if N % d == 0)
        s += chi3(N)*dN*N*mp.besselk(0, 2*mp.pi*N*n/3)
    return mp.mpf(4)/3*mp.power(n, 1.5)*s
for n in [1, 2, 4, 13, 41]:
    cf = phiphi(n); fv = r0['g0'] if n == 1 else r0['gv'][n]
    check('Lemma 2 closed form (own besselk) vs Fourier inversion of Xi^2 at log %d' % n, abs(cf - fv) < mp.mpf(10)**-50, '%s vs %s' % (mp.nstr(cf, 20), mp.nstr(fv, 20)))
check('B=1: identity residual < 1e-45 (Lemma 1, constants log 3, gamma + log pi, kernel)', abs(r0['res']) < mp.mpf(10)**-45, mp.nstr(r0['res'], 3))
check('B=1: A_1 converged (48 vs 96 GL nodes per panel)', r0['A1diff'] < mp.mpf(10)**-50, mp.nstr(r0['A1diff'], 2))
At, g0t = A1_tside(b0)
check('B=1: A_1 t-side (Gauss integral route) = A_1 u-side to 1e-25', abs(At - r0['A1']) < mp.mpf(10)**-25, '%s vs %s' % (mp.nstr(At, 20), mp.nstr(r0['A1'], 20)))
check('B=1: g(0) by mp.quad = trapezoid to 1e-25', abs(g0t - r0['g0']) < mp.mpf(10)**-25, mp.nstr(g0t, 20))
# aliasing check: the subgrid of step 2h must give the same g at u <= log 400 (2 pi/(2h) = 12.6)
gh2 = r0['gh']
def gfun2(gh, u):
    return 2*h/(2*mp.pi)*(gh[0] + 2*mp.fsum(gh[2*k]*mp.cos(2*k*h*u) for k in range(1, K//2 + 1)))
dmax = max(abs(gfun2(gh2, mp.log(n)) - r0['gv'][n]) for n in [2, 3, 5, 13, 41, 97])
check('trapezoid step h vs 2h agree at log 2 .. log 97 (aliasing negligible)', dmax < mp.mpf(10)**-50, mp.nstr(dmax, 2))
check('B=1: Defect/g(0) in [4.531e-2, 4.532e-2] (table)', mp.mpf('4.531e-2') < r0['ratio'] < mp.mpf('4.532e-2'), mp.nstr(r0['ratio'], 12))
check('B=1: 3 - bound = 0.2600 (table), bound 2.740045', abs(3 - r0['bound'] - mp.mpf('0.25995')) < 0.00001 and abs(r0['bound'] - mp.mpf('2.740045')) < 1e-6, mp.nstr(r0['bound'], 12))
check('B=1: largest share p = 3', r0['top'][0] == 3)

# ---------------------------------------------------------------- 2. the stored certificates
TABLE = {0: ('4.531e-2', '4.532e-2', '0.2600'), 2: ('6.130e-3', '6.131e-3', '3.656e-2'), 4: ('4.470e-4', '4.471e-4', '2.681e-3'),
         6: ('1.380e-5', '1.381e-5', '8.284e-5'), 8: ('5.136e-7', '5.137e-7', '3.082e-6'), 10: ('2.751e-9', '2.752e-9', '1.651e-8'),
         12: ('7.479e-11', '7.480e-11', '4.488e-10'), 14: ('5.252e-13', '5.253e-13', '3.152e-12'), 16: ('1.510e-17', '1.511e-17', '9.062e-17'),
         18: ('1.509e-17', '1.510e-17', '9.060e-17'), 20: ('1.509e-17', '1.510e-17', '9.060e-17'), 22: ('6.349e-23', '6.350e-23', '3.810e-22')}
CARRIER = {0: 3, 2: 5, 4: 7, 6: 3, 8: 13, 10: 17, 12: 19, 14: 23, 16: 31, 18: 31, 20: 31, 22: 41}
rows = {}
for deg in sorted(TABLE):
    b = load(deg)
    r = run(b, 'deg %d' % deg)
    rows[deg] = r
    lo, up, gap = TABLE[deg]
    check('deg %d: identity residual < 1e-45 g(0)' % deg, abs(r['res']) < mp.mpf(10)**-45*r['g0'], mp.nstr(r['res'], 3))
    check('deg %d: Defect/g(0) in [%s, %s]' % (deg, lo, up), mp.mpf(lo) < r['ratio'] < mp.mpf(up), mp.nstr(r['ratio'], 10))
    check('deg %d: 3 - bound rounds to %s' % (deg, gap), abs((3 - r['bound'])/mp.mpf(gap) - 1) < 2e-4, mp.nstr(3 - r['bound'], 10))
    check('deg %d: largest share at p = %d' % (deg, CARRIER[deg]), r['top'][0] == CARRIER[deg], str(r['top']))
    check('deg %d: two minimum methods agree to 1e-40' % deg, r['mindis'] < mp.mpf(10)**-40, mp.nstr(r['mindis'], 2))
    # B > 0 on R, independently of sympy: no root of q(s) = sum b_k s^k on [0, oo) (mpmath polyroots, 60 digits; computed)
    if len(b) > 1:
        rq = mp.polyroots([mp.mpf(x.numerator)/x.denominator for x in reversed(b)], maxsteps=500, extraprec=300)
        bad = [x for x in rq if abs(mp.im(x)) < mp.mpf(10)**-30*abs(x) and mp.re(x) >= 0]
        check('deg %d: q(s) has no root on [0, oo) (polyroots; B > 0 on R)' % deg, not bad and b[0] > 0,
              'roots nearest the positive axis: %s' % mp.nstr(min(abs(mp.arg(x)) for x in rq), 4))
    check('deg %d: p = 2 carries no cost (Defect_2/g(0) < 1e-50)' % deg, r['per'][2][0]/r['g0'] < mp.mpf(10)**-50, mp.nstr(r['per'][2][0]/r['g0'], 3))
    print('   g(log 2), g(log 4) / g(0) = %s, %s' % (mp.nstr(r['gv'][2]/r['g0'], 6), mp.nstr(r['gv'][4]/r['g0'], 6)))
print('deg 16: g(log n)/g(0), n = 2, 3, 4, 5, 7: %s' % ', '.join(mp.nstr(rows[16]['gv'][n]/rows[16]['g0'], 4) for n in [2, 3, 4, 5, 7]))
At8, g08 = A1_tside(load(8))
check('deg 8: A_1 t-side = u-side to 1e-25', abs(At8 - rows[8]['A1']) < mp.mpf(10)**-25, '%s vs %s' % (mp.nstr(At8, 20), mp.nstr(rows[8]['A1'], 20)))
check('deg 8: g(0) by mp.quad = trapezoid to 1e-25', abs(g08 - rows[8]['g0']) < mp.mpf(10)**-25, mp.nstr(g08, 20))

# ---------------------------------------------------------------- 3. Lemma 1 at B the author did not use (incl. B not >= 0)
for name, b in [('B=t^2', [Fraction(0), Fraction(1)]), ('B=1-t^2/10 (not >= 0)', [Fraction(1), Fraction(-1, 10)]),
                ('B=(t^2-64.64)^2', [Fraction(6464, 100)**2, Fraction(-2*6464, 100), Fraction(1)])]:
    r = run(b, name, full=False)
    check('%s: identity residual < 1e-45 |g|' % name, abs(r['res']) < mp.mpf(10)**-45*max(abs(r['g0']), 1e-10), mp.nstr(r['res'], 3))

# ---------------------------------------------------------------- 4. negative control: ghat = Xi_{-4}^2 in the chi_{-3} formula
XI4, w4 = xi_grid(4, (0, 1, 0, -1))
check('Xi_{-4} real on the grid', w4 < mp.mpf(10)**-45, mp.nstr(w4, 3))
r4 = run([Fraction(1)], 'control ghat = Xi_{-4}^2', XI=XI4, full=False)
# zeros of Xi_{-3} on (0, 70): sign changes on the grid, refined
zs = []
for k in range(K):
    if k*h > 70: break
    if XI3[k]*XI3[k + 1] < 0:
        with mp.workdps(45):
            zs.append(mp.findroot(lambda t: Lam(mp.mpf(1)/2 + 1j*t).real, (k*h, (k + 1)*h), solver='anderson'))
# smooth count: theta(T)/pi with theta = arg(3^{iT/2} pi^{-iT/2} Gamma(3/4 + iT/2)), continuous
def theta(T): return T/2*mp.log(3/mp.pi) + mp.im(mp.loggamma(mp.mpf(3)/4 + 1j*T/2)) - mp.im(mp.loggamma(mp.mpf(3)/4))
print('zeros of L(s, chi_-3) in (0, 70): %d (first %s); theta(70)/pi = %s' % (len(zs), mp.nstr(zs[0], 12), mp.nstr(theta(70)/mp.pi, 6)), flush=True)
check('zero count in (0,70) within 1.5 of theta(70)/pi (no missed pair)', abs(len(zs) - theta(70)/mp.pi) < 1.5)
Z = 2*mp.fsum(Lam(mp.mpf(1)/2 + 1j*z, 4, (0, 1, 0, -1)).real**2 for z in zs)
check('control: residual (log 3) g(0) + A_1 - 2 sum = zero sum 2 sum Xi_{-4}(gamma)^2 (to 1e-30)', abs(r4['res'] - Z) < mp.mpf(10)**-30,
      'residual %s, zero sum %s' % (mp.nstr(r4['res'], 20), mp.nstr(Z, 20)))
check('control: residual is far from 0 (the identity test is not vacuous)', abs(r4['res']) > mp.mpf(10)**-12, mp.nstr(r4['res'], 5))

print('\n%d checks, %d failed   [%.0f s]' % (NCHK[0], NCHK[1], time.time() - T0))
