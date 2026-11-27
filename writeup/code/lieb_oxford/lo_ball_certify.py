#!/usr/bin/env python3
"""Certified upper bound for the Lieb-Oxford constant with mu = nu = uniform unit ball
and the admissible function f = g of Lewin-Lieb-Seiringer (Lett. Math. Phys. 2022,
arXiv 2203.12473v3), Lemma 10, eq. (51).

    c_LO <= (3/2) (2 D(B,B)^2 J)^(1/3),   D(B,B) = 3/5,   J = 4 pi int_0^oo g(r) r^-5 dr.

Parts
 1. exact (sympy) proof of the closed formula (78) for Psi_BB and of the sphere formula (43)
 2. numerical cross-checks of (78) by two independent routes (Fourier integral, sphere mixture)
 3. exact algebra used below (region-A form, diagonal, factorisations, kappa)
 4. computer-assisted (arb ball arithmetic) proof that g(a) = Psi(a,a)/2 on [a0, a1]
 5. S0 (rigorous analytic bound on [0,a0]) and S1 (exact integral on [a0,a1])
 6. S2 = certified upper bound of int_{a1}^{R} g r^-5 by Taylor-in-a + convex chord + branch and bound in b
 7. S3 = rigorous tail bound on [R, oo)
 8. assembly, outward rounding, float cross-check
Run:  python3 lo_ball_certify.py [rho]      (rho = relative a-step in part 6, default 1e-3)
"""
import sys, time, math, heapq
import sympy as sp
from flint import arb, arb_poly, fmpq, ctx

ctx.prec = 128
T0 = time.time()
NCHK = [0, 0]


def check(cond, msg):
    NCHK[0] += 1
    if cond:
        print('ok   ', msg, flush=True)
    else:
        NCHK[1] += 1
        print('FAIL ', msg, flush=True)


a, b, u, s = sp.symbols('a b u s', positive=True)
R1, R2 = sp.symbols('R1 R2', positive=True)

# ----------------------------------------------------------------------------------------------
# Part 1. Exact proof of (78) (ball) and (43) (sphere), d = 1, R1 = 1/a, R2 = 1/b.
# The interaction W = 2 D(mu_{0,a}, nu_{e1,b}) of two unit charges is computed as
#     W = int_0^oo p(s) V1(s) ds,
# V1 = potential of the charge at 0 (radius R1), p = density of |y| for y distributed as the
# charge centred at e1 (radius R2).  Phi = a^3 b^3 (1 - W).
# ----------------------------------------------------------------------------------------------
print('--- Part 1: exact proof of (78) and (43)')
Vin = (3*R1**2 - s**2)/(2*R1**3)            # uniform ball, |x| = s <= R1
Vout = 1/s                                  # |x| >= R1 (Newton)
pin = 3*s**2/R2**3                          # |y| = s, sphere {|y|=s} inside B(e1,R2): s <= R2 - 1
pcap = 3*s*(R2**2 - (s - 1)**2)/(4*R2**3)   # |1-R2| <= s <= 1+R2 (spherical cap area / volume)
check(sp.simplify(sp.integrate(pcap, (s, 1 - R2, 1 + R2)) - 1) == 0, 'ball density normalised (R2<=1)')
check(sp.simplify(sp.integrate(pin, (s, 0, R2 - 1)) + sp.integrate(pcap, (s, R2 - 1, R2 + 1)) - 1) == 0,
      'ball density normalised (R2>=1)')


def T1f(x, y):
    return (x + y - x*y)**4/160*(x**2*y**2 - 5*x**2 - 5*y**2 + 4*x*y**2 + 4*y*x**2 + 20*x*y)


def T2f(x, y, amb):  # amb = |x - y|
    return (amb - x*y)**4/160*(x**2*y**2 - 5*x**2 - 5*y**2 + 4*amb*x*y - 20*x*y)


A_, B_ = 1/R1, 1/R2
PHI = lambda W: (1 - W)/(R1**3*R2**3)
cases = []
# (I) disjoint, R1 + R2 <= 1: (78) gives 0
W = sp.integrate(pcap*Vout, (s, 1 - R2, 1 + R2))
cases.append(('(I) disjoint', sp.simplify(PHI(W))))
# (III) partial overlap |R1-R2| < 1 < R1+R2, R2 <= 1 and R2 >= 1: (78) = T1 only
W = sp.integrate(pcap*Vin, (s, 1 - R2, R1)) + sp.integrate(pcap*Vout, (s, R1, 1 + R2))
cases.append(('(III) overlap, R2<=1', sp.simplify(sp.expand(PHI(W) - T1f(A_, B_)))))
W = (sp.integrate(pin*Vin, (s, 0, R2 - 1)) + sp.integrate(pcap*Vin, (s, R2 - 1, R1))
     + sp.integrate(pcap*Vout, (s, R1, 1 + R2)))
cases.append(('(III) overlap, R2>=1', sp.simplify(sp.expand(PHI(W) - T1f(A_, B_)))))
# (II) ball 1 inside ball 2, R2 >= R1 + 1 (a > b)
W_II = (sp.integrate(pin*Vin, (s, 0, R1)) + sp.integrate(pin*Vout, (s, R1, R2 - 1))
        + sp.integrate(pcap*Vout, (s, R2 - 1, R2 + 1)))
cases.append(('(II) 1 inside 2', sp.simplify(sp.expand(PHI(W_II) - (T1f(A_, B_) - T2f(A_, B_, A_ - B_))))))
cases.append(('(II) closed W', sp.simplify(W_II - (3*R2**2 - 1 - sp.Rational(3, 5)*R1**2)/(2*R2**3))))
# (II') ball 2 inside ball 1, R1 >= R2 + 1 (b > a)
W = sp.integrate(pcap*Vin, (s, 1 - R2, 1 + R2))
cases.append(("(II') 2 inside 1, R2<=1", sp.simplify(sp.expand(PHI(W) - (T1f(A_, B_) - T2f(A_, B_, B_ - A_))))))
W = sp.integrate(pin*Vin, (s, 0, R2 - 1)) + sp.integrate(pcap*Vin, (s, R2 - 1, R2 + 1))
cases.append(("(II') 2 inside 1, R2>=1", sp.simplify(sp.expand(PHI(W) - (T1f(A_, B_) - T2f(A_, B_, B_ - A_))))))
for name, d in cases:
    check(d == 0, 'formula (78) in configuration ' + name)
# in each configuration the inactive (.)_+ factors vanish: (I) x1 = ab(R1+R2-1) <= 0; (III) x2 = ab(|R1-R2|-1) < 0;
# (II),(II') x2 >= 0 and x1 >= x2.  These sign facts are elementary and are restated in the note.

# sphere formula (43): p_sigma(s) = s/(2 R2) on [|1-R2|, 1+R2], V1 = 1/max(s, R1)
ps = s/(2*R2)
SIG = lambda x, y, amb: x**2*y**2/4*(sp.Max(x + y - x*y, 0)**2 - sp.Max(amb - x*y, 0)**2)
sig_cases = []
W = sp.integrate(ps/s, (s, 1 - R2, 1 + R2))       # R2 <= 1, R1 <= 1-R2 : disjoint
sig_cases.append(('disjoint', sp.simplify(PHI(W))))
W = sp.integrate(ps/s, (s, R2 - 1, 1 + R2))       # R2 >= 1, R1 <= R2-1 : sphere 1 inside sphere 2
sig_cases.append(('1 inside 2', sp.simplify(sp.expand(PHI(W) - A_**2*B_**2/4*((A_ + B_ - A_*B_)**2 - (A_ - B_ - A_*B_)**2)))))
W = sp.integrate(ps/R1, (s, 1 - R2, 1 + R2))      # R1 >= 1+R2, R2 <= 1 : sphere 2 inside sphere 1
sig_cases.append(('2 inside 1 (R2<=1)', sp.simplify(sp.expand(PHI(W) - A_**2*B_**2/4*((A_ + B_ - A_*B_)**2 - (B_ - A_ - A_*B_)**2)))))
W = sp.integrate(ps/R1, (s, R2 - 1, 1 + R2))      # R1 >= 1+R2, R2 >= 1
sig_cases.append(('2 inside 1 (R2>=1)', sp.simplify(sp.expand(PHI(W) - A_**2*B_**2/4*((A_ + B_ - A_*B_)**2 - (B_ - A_ - A_*B_)**2)))))
W = sp.integrate(ps/R1, (s, 1 - R2, R1)) + sp.integrate(ps/s, (s, R1, 1 + R2))   # overlap, R2 <= 1
sig_cases.append(('overlap (R2<=1)', sp.simplify(sp.expand(PHI(W) - A_**2*B_**2/4*(A_ + B_ - A_*B_)**2))))
W = sp.integrate(ps/R1, (s, R2 - 1, R1)) + sp.integrate(ps/s, (s, R1, 1 + R2))   # overlap, R2 >= 1
sig_cases.append(('overlap (R2>=1)', sp.simplify(sp.expand(PHI(W) - A_**2*B_**2/4*(A_ + B_ - A_*B_)**2))))
for name, d in sig_cases:
    check(d == 0, 'sphere formula (43) in configuration ' + name)

# ----------------------------------------------------------------------------------------------
# Part 2. Numerical cross-checks of (78) by two independent routes.
# ----------------------------------------------------------------------------------------------
print('--- Part 2: independent numerical routes for (78)')
import mpmath as mp
import numpy as np
from scipy.integrate import dblquad


def psi78_float(x, y):
    x1 = x + y - x*y
    x2 = abs(x - y) - x*y
    v = 0.0
    if x1 > 0:
        v += x1**4/160*(x*x*y*y - 5*x*x - 5*y*y + 4*x*y*y + 4*y*x*x + 20*x*y)
    if x2 > 0:
        v -= x2**4/160*(x*x*y*y - 5*x*x - 5*y*y + 4*abs(x - y)*x*y - 20*x*y)
    return v


def psi_fourier(x, y):
    # W = (2/pi) int_0^oo F(k R1) F(k R2) sin(k)/k dk, F(z) = 3 (sin z - z cos z)/z^3 (unit-charge ball)
    mp.mp.dps = 25
    r1, r2 = mp.mpf(1)/x, mp.mpf(1)/y
    F = lambda z: 3*(mp.sin(z) - z*mp.cos(z))/z**3 if abs(z) > mp.mpf('1e-4') else 1 - z**2/10 + z**4/280
    f = lambda k: F(k*r1)*F(k*r2)*mp.sin(k)/k
    W = 2/mp.pi*mp.quadosc(f, [0, mp.inf], omega=1)
    return float(mp.mpf(x)**3*mp.mpf(y)**3*(1 - W))


def psi_spheres(x, y):
    # uniform ball = int_0^1 3 r^2 dr sigma_r ; Phi_{sigma_r sigma_t}(x,y) = r^3 t^3 Phi_ss(x/r, y/t) with (43)
    def pss(p, q):
        return p*p*q*q/4*(max(p + q - p*q, 0)**2 - max(abs(p - q) - p*q, 0)**2)
    f = lambda t, r: pss(x/r, y/t)*r**3*t**3*9*r*r*t*t
    v, e = dblquad(f, 0, 1, 0, 1, epsabs=1e-13, epsrel=1e-11)
    return v


test_pts = [(0.3, 0.2), (1.0, 0.7), (2.5, 0.6), (1.7, 1.2), (0.8, 0.8), (5.0, 1.1), (4.0, 0.3), (0.05, 0.01),
            (0.889, 0.889), (3.0, 1.5 + 1e-3), (3.0, 1.5 - 1e-3),          # next to a+b-ab = 0
            (3.0, 0.75 + 1e-3), (3.0, 0.75 - 1e-3),                        # next to a-b-ab = 0
            (1.0, 1.0 - 1e-6), (30.0, 0.618034), (20.0, 1.02)]
for (x, y) in test_pts:
    v78, vF = psi78_float(x, y), psi_fourier(x, y)
    check(abs(v78 - vF) <= 1e-12*max(1.0, x**3*y**3), f'(78) = Fourier route at (a,b)=({x},{y}): {v78:.15e} vs {vF:.15e}')
for (x, y) in [(1.0, 0.7), (2.5, 0.6), (5.0, 1.1), (0.8, 0.8)]:
    v78, vS = psi78_float(x, y), psi_spheres(x, y)
    check(abs(v78 - vS) <= 1e-8*max(1.0, x**3*y**3), f'(78) = sphere-mixture route at ({x},{y}): {v78:.12e} vs {vS:.12e}')

# ----------------------------------------------------------------------------------------------
# Part 3. Exact algebra used below.  For b <= a:
#   x1 = a + b - ab, x2 = a - b - ab (x2 <= x1);  Psi = PA if x2 >= 0, T1 if x2 <= 0 <= x1, 0 if x1 <= 0,
#   PA = T1 - T2 = a^3 b^3 (1 - 3b/2 + b^3/2) + 3 a b^6/10   (region A: ball 1/a inside ball 1/b)
#   Delta(b) = Psi(b,b)/2 = b^6 (2-b)_+^4 (b^2 + 8b + 10)/320,   h(a,b) = Psi(a,b) - Delta(b).
# ----------------------------------------------------------------------------------------------
print('--- Part 3: exact algebra')
T1 = sp.expand(T1f(a, b))
T2 = sp.expand(T2f(a, b, a - b))
PA = sp.expand(a**3*b**3*(1 - sp.Rational(3, 2)*b + b**3/2) + sp.Rational(3, 10)*a*b**6)
DL = sp.expand(b**6*(2 - b)**4*(b**2 + 8*b + 10)/320)
check(sp.expand(T1 - T2 - PA) == 0, 'T1 - T2 = PA identically')
check(sp.expand(T1.subs(a, b)/2 - DL) == 0, 'Psi(b,b)/2 = Delta(b) for b <= 2 (x2 = -b^2 < 0)')
# region A with a -> oo: Psi/a^3 -> k(b) = b^3 (1 - 3b/2 + b^3/2); max on [0,1] at the golden ratio point
kb = b**3*(1 - sp.Rational(3, 2)*b + b**3/2)
check(sp.factor(sp.diff(kb, b)) == sp.factor(3*b**2*(b - 1)*(b**2 + b - 1)), "k'(b) = 3 b^2 (b-1)(b^2+b-1)")
bstar = (sp.sqrt(5) - 1)/2
kappa_exact = sp.nsimplify(sp.simplify(kb.subs(b, bstar)))
check(sp.simplify(kappa_exact - (5*sp.sqrt(5) - 11)/4) == 0, f'kappa = max k = (5 sqrt5 - 11)/4 = {float(kappa_exact):.15f}')
# triangle quantities, b = a u
DT_B = sp.expand(DL.subs(b, a) - T1.subs(b, a*u) + DL.subs(b, a*u))
DT_A = sp.expand(DL.subs(b, a) - PA.subs(b, a*u) + DL.subs(b, a*u))
qB, rB = sp.div(sp.Poly(DT_B, u), sp.Poly(sp.expand((1 - u)**2*a**6), u))
check(rB.is_zero, 'D_B(a,u) = a^6 (1-u)^2 q_B(a,u) exactly')
qB = sp.expand(qB.as_expr())
EA = sp.expand(DT_A/a**6)
qA, rA = sp.div(sp.Poly(EA, u), sp.Poly(sp.expand((1 - u)**2), u))
qA, rA = sp.expand(qA.as_expr()), sp.expand(rA.as_expr())
check(sp.expand(qA*(1 - u)**2 + rA - EA) == 0, 'D_A/a^6 = (1-u)^2 q_A + r_A')
qB1 = sp.factor(qB.subs(u, 1))
check(sp.expand(qB1 - sp.Rational(9, 80)*(a - 2)**2*(a**4 + 4*a**3 - 8*a**2 - 8*a + 10)) == 0,
      'q_B(a,1) = 9 (a-2)^2 (a^4+4a^3-8a^2-8a+10)/80')
hb = T1 - DL
check(sp.expand(sp.diff(hb, b).subs(b, a)) == 0, 'd/db h(a,b) = 0 on the diagonal b = a')
check(sp.expand(sp.diff(hb, b, 2).subs(b, a) + 2*a**4*qB1) == 0, 'h_bb(a,a) = -2 a^4 q_B(a,1)')
ac_poly = sp.Poly(a**4 + 4*a**3 - 8*a**2 - 8*a + 10, a)
ac_roots = sorted([float(sp.re(r)) for r in ac_poly.nroots(n=30) if abs(sp.im(r)) < 1e-20 and sp.re(r) > 0])
print('    positive real roots of a^4+4a^3-8a^2-8a+10:', ac_roots)

# ---------- conversion of sympy polynomials to arb machinery ----------


def to_fmpq(c):
    c = sp.Rational(c)
    return fmpq(int(c.p), int(c.q))


class BiPoly:
    """P(a,b) with exact rational coefficients; evaluate as arb_poly in b at an arb value of a."""

    def __init__(self, expr, xa=a, xb=b):
        P = sp.Poly(sp.expand(expr), xb, xa)
        self.degb = P.degree(xb)
        self.dega = P.degree(xa)
        self.C = [[fmpq(0)]*(self.dega + 1) for _ in range(self.degb + 1)]
        for (j, i), c in P.terms():
            self.C[j][i] = to_fmpq(c)

    def in_b(self, A):
        powA = [arb(1)]
        for _ in range(self.dega):
            powA.append(powA[-1]*A)
        return arb_poly([sum((powA[i]*self.C[j][i] for i in range(self.dega + 1) if self.C[j][i] != 0), arb(0))
                         for j in range(self.degb + 1)])


def _hull2(vals):
    x = vals[0]
    for v in vals[1:]:
        x = x.union(v)
    return x


# ----------------------------------------------------------------------------------------------
# Part 4. g(a) = Psi(a,a)/2 on [a0, a1]: D(a,u) = Delta(a) - Psi(a,au) + Delta(au) >= 0 for u in [0,1].
# Region A: u <= 1/(1+a), D/a^6 = (1-u)^2 q_A + r_A.  Region B: u >= 1/(1+a), D/a^6 = (1-u)^2 q_B.
# (For a <= 1, x1 = a(1 + u(1-a)) > 0, so no zero region.)
# ----------------------------------------------------------------------------------------------
print('--- Part 4: g = Psi(a,a)/2 on [a0, a1] (ball arithmetic)')
A0_ = fmpq(1, 1024)
A1_ = fmpq(7, 8)         # 0.875 < a_c = 0.88904...
check(float(A1_) < ac_roots[0], f'a1 = 7/8 < a_c = {ac_roots[0]:.12f}')
qA_bp, rA_bp, qB_bp = BiPoly(qA, a, u), BiPoly(rA, a, u), BiPoly(qB, a, u)


def tri_box_ok(alo, ahi, ulo, uhi):
    Abox = arb(alo).union(arb(ahi))
    ok = True
    # region A part: u in [ulo, min(uhi, 1/(1+alo))]
    ucut = (1/(1 + arb(alo))).upper()
    if arb(ulo) <= ucut:
        uh = arb(uhi) if arb(uhi) < ucut else ucut
        U = arb(ulo).union(arb(uh))
        val = (1 - U)**2*qA_bp.in_b(Abox)(U) + rA_bp.in_b(Abox)(U)
        ok = ok and bool(val > 0)
    # region B part: u in [max(ulo, 1/(1+ahi)), uhi]
    ucut2 = (1/(1 + arb(ahi))).lower()
    if arb(uhi) >= ucut2:
        ul = arb(ulo) if arb(ulo) > ucut2 else ucut2
        U = ul.union(arb(uhi))
        val = qB_bp.in_b(Abox)(U)
        ok = ok and bool(val > 0)
    return ok


def tri_verify(alo, ahi, ulo, uhi, depth=0, stats=[0, 0]):
    stats[0] += 1
    if tri_box_ok(alo, ahi, ulo, uhi):
        return True
    if depth > 40:
        stats[1] += 1
        return False
    am, um = (alo + ahi)/2, (ulo + uhi)/2
    return all(tri_verify(x0, x1, y0, y1, depth + 1, stats) for (x0, x1) in ((alo, am), (am, ahi))
               for (y0, y1) in ((ulo, um), (um, uhi)))


st = [0, 0]
okT = True
NA, NU = 64, 16
for i in range(NA):
    alo = A0_ + (A1_ - A0_)*fmpq(i, NA)
    ahi = A0_ + (A1_ - A0_)*fmpq(i + 1, NA)
    for j in range(NU):
        okT = okT and tri_verify(alo, ahi, fmpq(j, NU), fmpq(j + 1, NU), 0, st)
check(okT, f'D(a,u) > 0 on [1/1024, 7/8] x [0,1] except the zero line u = 1 ({st[0]} boxes)')

# ----------------------------------------------------------------------------------------------
# Part 5. S0 and S1.
# S0: for a <= a0: Psi(a,b) <= a^3 b^3 (W >= 0), Psi(b,b) >= b^6 (1 - 3b/2) (eq. (21): W <= b V_B(b e1)
#     = b(3-b^2)/2), and a^3b^3 - b^6(1-3b/2)/2 is increasing in b on [0,a]  =>  g(a) <= a^6/2 + 3a^7/4.
#     int_0^a0 (r/2 + 3 r^2/4) dr = a0^2/4 + a0^3/4.
# S1: int_a0^a1 Psi(r,r)/(2 r^5) dr exact.
# ----------------------------------------------------------------------------------------------
print('--- Part 5: S0, S1')
check(sp.expand(sp.diff(a**3*b**3 - b**6*(1 - sp.Rational(3, 2)*b)/2, b) - 3*b**2*(a**3 - b**3)
                - sp.Rational(21, 4)*b**6) == 0, 'd/db [a^3b^3 - b^6(1-3b/2)/2] = 3b^2(a^3-b^3) + 21 b^6/4 >= 0')
a0r, a1r = sp.Rational(1, 1024), sp.Rational(7, 8)
S0 = a0r**2/4 + a0r**3/4
S1 = sp.integrate(sp.expand(DL.subs(b, a)/a**5), (a, a0r, a1r))
check(S1.is_rational, f'S1 exact rational = {float(S1):.15e}')
S0a, S1a = arb(to_fmpq(S0)), arb(to_fmpq(S1))

# ----------------------------------------------------------------------------------------------
# Part 6. S2 = int_{a1}^{R} g(r) r^-5 dr.
# For a in I = [A0,A1], am = mid, d = A1 - am, and every b in [0, A1]:
#   H(a,b) <= H(am,b) + (a-am) H_a(am,b) + (1/2) Haa+ (a-am)^2,   Haa+ >= sup_{I x [0,A1]} H_aa, Haa+ >= 0.
# G(s) = max_{b in [0,A1]} [H(am,b) + s H_a(am,b)] is convex in s, so G(s) <= chord through G(-d), G(d).
# Hence g(a) <= alpha + beta (a-am) + Haa+ d^2/2 on I with alpha, beta from certified upper bounds of G(+-d).
# G(+-d): branch and bound in b with the Taylor box bound F(c) + |F'(c)| r + max(F''(box),0) r^2/2.
# ----------------------------------------------------------------------------------------------
print('--- Part 6: S2 by Taylor/chord/branch-and-bound (ball arithmetic)')
RHO = float(sys.argv[1]) if len(sys.argv) > 1 else 1e-3
RTAIL = 100
PA_bp, T1_bp = BiPoly(PA), BiPoly(T1)
PAa_bp, T1a_bp = BiPoly(sp.diff(PA, a)), BiPoly(sp.diff(T1, a))
PAaa_bp, T1aa_bp = BiPoly(sp.diff(PA, a, 2)), BiPoly(sp.diff(T1, a, 2))
DLp = arb_poly([arb(to_fmpq(c)) for c in reversed(sp.Poly(DL, b).all_coeffs())])
DLp1, DLp2 = DLp.derivative(), DLp.derivative().derivative()
TWO = arb(2)


def pieces(Aval, Bval):
    """which of the pieces A (x2>=0), B (x2<=0<=x1), Z (x1<=0) can occur for a in Aval, b in Bval"""
    x1 = Aval + Bval - Aval*Bval
    x2 = Aval - Bval - Aval*Bval
    out = []
    if not bool(x2 < 0):
        out.append('A')
    if (not bool(x2 > 0)) and (not bool(x1 < 0)):
        out.append('B')
    if not bool(x1 > 0):
        out.append('Z')
    return out


def dl_part(Bval, k):
    p = (DLp, DLp1, DLp2)[k]
    if bool(Bval < TWO):
        return p(Bval)
    if bool(Bval > TWO):
        return arb(0)
    return _hull2([p(Bval), arb(0)])


class Slice:
    """F(b) = H(am,b) + sv H_a(am,b) and its b-derivatives, piecewise."""

    def __init__(self, am, sv):
        self.am = am
        FA = PA_bp.in_b(am) + PAa_bp.in_b(am)*sv
        FB = T1_bp.in_b(am) + T1a_bp.in_b(am)*sv
        self.P = {'A': [FA, FA.derivative(), FA.derivative().derivative()],
                  'B': [FB, FB.derivative(), FB.derivative().derivative()]}

    def ev(self, Bval, k):
        vals = []
        for pc in pieces(self.am, Bval):
            vals.append(arb(0) if pc == 'Z' else self.P[pc][k](Bval))
        psi = _hull2(vals)
        return psi - dl_part(Bval, k)


def bb_max(sl, bmax, tol_rel=1e-10, maxit=20000):
    """certified upper bound of max_{b in [0,bmax]} F(b); bmax is a float (exact dyadic)"""
    am = float(sl.am.mid())
    cuts = sorted(set([0.0, bmax] + [bmax*k/16 for k in range(1, 16)] +
                      [x for x in (am/(1 + am), am/(am - 1) if am > 1 else -1, 2.0) if 0 < x < bmax]))
    L = [None]
    heap = []

    def push(lo, hi):
        cb = (arb(lo) + arb(hi))/2      # exact centre of [lo,hi] (128-bit arb, lo/hi are doubles)
        r = (arb(hi) - arb(lo))/2       # exact half width
        Fc = sl.ev(cb, 0)
        F1 = sl.ev(cb, 1)
        F2 = sl.ev(arb(lo).union(arb(hi)), 2)
        F2u = F2.upper()
        U = (Fc.upper() + abs(F1).upper()*r + (F2u if F2u > 0 else arb(0))*(r*r/2)).upper()
        lowc = Fc.lower()
        if L[0] is None or lowc > L[0]:
            L[0] = lowc
        heapq.heappush(heap, (-float(U), lo, hi, U))
    for lo, hi in zip(cuts[:-1], cuts[1:]):
        push(lo, hi)
    it = 0
    while True:
        negU, lo, hi, U = heap[0]
        tol = tol_rel*max(abs(float(L[0])), 1e-30)
        if float(U) <= float(L[0]) + tol or it >= maxit or hi - lo < 1e-15*max(1.0, bmax):
            best = heap[0][3]
            for e in heap:               # exact maximum of all box bounds (robust to float ordering)
                if e[3] > best:
                    best = e[3]
            return best, it, float(L[0])
        heapq.heappop(heap)
        m = (lo + hi)/2
        push(lo, m)
        push(m, hi)
        it += 1


# factored second a-derivatives (checked against sympy):
#   PA_aa = 6 a b^3 (1 - 3b/2 + b^3/2)
#   T1_aa = [12 x1^2 (1-b)^2 P1 + 8 x1^3 (1-b) P1_a + x1^4 P1_aa]/160,  x1 = a + b - ab
P1s = a**2*b**2 - 5*a**2 - 5*b**2 + 4*a*b**2 + 4*a**2*b + 20*a*b
x1s = a + b - a*b
T1aa_f = (12*x1s**2*(1 - b)**2*P1s + 8*x1s**3*(1 - b)*sp.diff(P1s, a) + x1s**4*sp.diff(P1s, a, 2))/160
check(sp.expand(T1aa_f - sp.diff(T1, a, 2)) == 0, 'factored form of d^2 T1/da^2')
check(sp.expand(6*a*b**3*(1 - sp.Rational(3, 2)*b + b**3/2) - sp.diff(PA, a, 2)) == 0, 'factored form of d^2 PA/da^2')


def haa_box(A, B):
    vals = []
    for pc in pieces(A, B):
        if pc == 'Z':
            vals.append(arb(0))
        elif pc == 'A':
            vals.append(6*A*B**3*(1 - 3*B/2 + B**3/2))
        else:
            x1 = A + B*(1 - A)
            P1 = A*A*B*B - 5*A*A - 5*B*B + 4*A*B*B + 4*A*A*B + 20*A*B
            P1a = 2*A*B*B - 10*A + 4*B*B + 8*A*B + 20*B
            P1aa = 2*B*B + 8*B - 10
            vals.append((12*x1**2*(1 - B)**2*P1 + 8*x1**3*(1 - B)*P1a + x1**4*P1aa)/160)
    return _hull2(vals)


def haa_sup(P, Q, nb=1024, na=4):
    """upper bound of max(0, sup H_aa) over [P,Q] x [0,Q]  (H_aa = Psi_aa); P,Q doubles.
    For a >= P > 1 and b >= P/(P-1): x1 = a - b(a-1) <= P - b(P-1) <= 0, so Psi = 0 there (Newton)."""
    best = arb(0)
    bsup = Q
    if P > 1:
        bs = (arb(P)/(arb(P) - 1)).upper()
        if bs < arb(Q):
            bsup = float(bs) if arb(float(bs)) >= bs else float(bs) + abs(float(bs))*2**-50
    ae = [P] + [P + (Q - P)*i/na for i in range(1, na)] + [Q]
    be_ = [0.0] + [bsup*j/nb for j in range(1, nb)] + [bsup]
    for i in range(na):
        Aval = arb(ae[i]).union(arb(ae[i + 1]))
        for j in range(nb):
            Bval = arb(be_[j]).union(arb(be_[j + 1]))
            v = haa_box(Aval, Bval).upper()
            if v > best:
                best = v
    return best


# grid of a-intervals: geometric with ratio (1+RHO); endpoints are doubles, consecutive intervals share them
grid = [0.875]
while grid[-1] < RTAIL:
    grid.append(min(grid[-1]*(1 + RHO), float(RTAIL)))
BLK = 64
S2 = arb(0)
nit_tot = 0
t6 = time.time()
gap_max = 0.0
haa_cache = {}
for k in range(len(grid) - 1):
    A0, A1 = grid[k], grid[k + 1]
    blk = k//BLK
    if blk not in haa_cache:
        P = grid[blk*BLK]
        Q = grid[min((blk + 1)*BLK, len(grid) - 1)]
        haa_cache[blk] = haa_sup(P, Q)
    Haa = haa_cache[blk]
    am = (arb(A0) + arb(A1))/2      # exact midpoint
    d = arb(A1) - am                # exact half width
    Gm, it1, Lm = bb_max(Slice(am, -d), A1)
    Gp, it2, Lp = bb_max(Slice(am, d), A1)
    nit_tot += it1 + it2
    gap_max = max(gap_max, (float(Gm) - Lm)/max(abs(Lm), 1e-300), (float(Gp) - Lp)/max(abs(Lp), 1e-300))
    alpha = (Gm + Gp)/2
    beta = (Gp - Gm)/(2*d)
    A0b, A1b = arb(A0), arb(A1)
    i5 = (A0b**-4 - A1b**-4)/4          # int_I r^-5
    i4 = (A0b**-3 - A1b**-3)/3          # int_I r^-4
    S2 += alpha*i5 + beta*(i4 - am*i5) + Haa*d*d/2*i5
    if k % 2000 == 0:
        print(f'    a = {A0:.6f}  G-,G+ = {float(Gm):.10e}, {float(Gp):.10e}  Haa+ = {float(Haa):.4e}  '
              f'iters {it1},{it2}  S2 so far {float(S2.upper()):.12f}  t = {time.time()-t6:.1f}s', flush=True)
print(f'    {len(grid)-1} a-intervals, {nit_tot} bisections, max relative B&B gap {gap_max:.2e}, {time.time()-t6:.1f}s')
# sanity checks of the enclosures by random float sampling (not part of the proof)
import random
random.seed(12345)
f_PA, f_T1 = sp.lambdify((a, b), PA), sp.lambdify((a, b), T1)
f_PAa, f_T1a = sp.lambdify((a, b), sp.diff(PA, a)), sp.lambdify((a, b), sp.diff(T1, a))
f_PAaa, f_T1aa = sp.lambdify((a, b), sp.diff(PA, a, 2)), sp.lambdify((a, b), sp.diff(T1, a, 2))
f_DL = sp.lambdify(b, DL)


def piecef(fA, fB, x, y):
    x1, x2 = x + y - x*y, x - y - x*y
    return fA(x, y) if x2 >= 0 else (fB(x, y) if x1 >= 0 else 0.0)


bad = 0
for blk, Hv in haa_cache.items():
    P = grid[blk*BLK]; Qv = grid[min((blk + 1)*BLK, len(grid) - 1)]
    for _ in range(40):
        x, y = random.uniform(P, Qv), random.uniform(0, Qv)
        if piecef(f_PAaa, f_T1aa, x, y) > float(Hv)*(1 + 1e-9) + 1e-12:
            bad += 1
check(bad == 0, f'sampling: Psi_aa <= Haa+ at {40*len(haa_cache)} random points of the blocks')
bad = 0
for k in range(0, len(grid) - 1, 997):
    A0, A1 = grid[k], grid[k + 1]
    am = (A0 + A1)/2; dd = A1 - am
    Gm, _, _ = bb_max(Slice(arb(am), -arb(dd)), A1)
    Gp, _, _ = bb_max(Slice(arb(am), arb(dd)), A1)
    for sgn, Gv in ((-1, Gm), (1, Gp)):
        bs = [A1*j/4000 for j in range(4001)]
        fm = max(piecef(f_PA, f_T1, am, y) + sgn*dd*piecef(f_PAa, f_T1a, am, y) - (f_DL(y) if y < 2 else 0.0) for y in bs)
        if fm > float(Gv) + 1e-12*abs(fm):
            bad += 1
check(bad == 0, 'sampling: certified G(+-d) >= float maximum on a 4001-point b-grid (every 997th interval)')
check(grid[0] == 0.875 and grid[-1] == float(RTAIL) and all(x < y for x, y in zip(grid[:-1], grid[1:])),
      'a-grid is increasing and covers [7/8, R] exactly')

# ----------------------------------------------------------------------------------------------
# Part 7. Tail a >= R (t = 1/a <= 1/R): with beta0 = 3/4,
#   region A (b <= 1/(1+t)): Psi/a^3 = k(b) + (3/10) t^2 b^6 <= max(kappa + 0.3 t^2 beta0^6, k(beta0) + 0.3 t^2)
#   region B (1/(1+t) <= b <= 1/(1-t)): W >= 1/(1+t)  =>  Psi/a^3 <= t/((1-t)^3 (1+t))
#   region C (b >= 1/(1-t)): Psi = 0;  and h <= Psi since Psi(b,b) >= 0.
# If k(beta0) + 0.3/R^2 <= kappa (and the region B bound <= kappa) at a = R, then for all a >= R
#   g(a) <= kappa a^3 + 0.3 beta0^6 a,  and  int_R^oo g r^-5 <= kappa/R + 0.1 beta0^6/R^3.
# ----------------------------------------------------------------------------------------------
print('--- Part 7: tail')
kap = (5*arb(5).sqrt() - 11)/4
be = arb(fmpq(3, 4))
kbe = be**3*(1 - 3*be/2 + be**3/2)
tR = 1/arb(RTAIL)
check(bool(kbe + arb(fmpq(3, 10))*tR**2 < kap), f'k(3/4) + 0.3/R^2 = {float(kbe.mid()+0.3/RTAIL**2):.6f} < kappa')
check(bool(tR/((1 - tR)**3*(1 + tR)) < kap), f'region-B bound at t=1/R: {float((tR/((1 - tR)**3*(1 + tR))).mid()):.6f} < kappa')
check(bool(3*(kap - kbe)*arb(RTAIL)**2 > arb(fmpq(3, 10))*(1 - be**6)),
      'kappa a^3 + 0.3 b0^6 a - (k(b0) a^3 + 0.3 a) increasing for a >= R')
S3 = kap/RTAIL + arb(fmpq(1, 10))*be**6/arb(RTAIL)**3
print(f'    S3 <= {float(S3.upper()):.12e}')

# ----------------------------------------------------------------------------------------------
# Part 8. Assembly.
# ----------------------------------------------------------------------------------------------
print('--- Part 8: assembly')
# D(B,B) = (1/2) int V_B dB with V_B(s) = (3 - s^2)/2 on the unit ball
check(sp.integrate(3*s**2*(3 - s**2)/2, (s, 0, 1))/2 == sp.Rational(3, 5), 'D(B,B) = 3/5')
Jup = 4*arb.pi()*(S0a + S1a + S2 + S3)
D = arb(fmpq(3, 5))
c = arb(fmpq(3, 2))*(2*Jup*D*D).root(3)
Ju = float(Jup.upper())
cu = float(c.upper())
print(f'    S0 = {float(S0):.6e}, S1 = {float(S1):.12f}, S2 <= {float(S2.upper()):.12f}, S3 <= {float(S3.upper()):.12f}')
print(f'    J <= {Jup.upper()}')
print(f'    c_LO <= {c.upper()}')
c_round = math.ceil(cu*1e6)/1e6
check(c_round >= cu, f'outward rounding: c_LO <= {c_round:.6f}')
check(c_round < 1.64, 'certified value below 1.64')

# float cross-check of J (independent quadrature of a float g)
from scipy.optimize import minimize_scalar
from scipy.integrate import quad


def gfl(x, nb=3001):
    bs = np.linspace(0, x, nb)
    vals = [psi78_float(x, y) - 0.5*psi78_float(y, y) for y in bs]
    i = int(np.argmax(vals))
    r = minimize_scalar(lambda y: -(psi78_float(x, y) - 0.5*psi78_float(y, y)),
                        bounds=(bs[max(i - 1, 0)], bs[min(i + 1, nb - 1)]), method='bounded', options={'xatol': 1e-14})
    return max(-r.fun, vals[i])


pts = [0, 0.5, 0.889042374847991, 1.0, 1.5, 2.0, 3.0, 5, 10, 20, 50, 100]
Jf = sum(quad(lambda r: gfl(r)/r**5, lo, hi, limit=200, epsabs=1e-13, epsrel=1e-11)[0] for lo, hi in zip(pts[:-1], pts[1:]))
Jf = 4*math.pi*(Jf + float(kappa_exact)/100)
cf = 1.5*(2*Jf*0.36)**(1/3)
print(f'    float cross-check: J = {Jf:.10f}, c = {cf:.10f}')
check(Jf <= Ju and Ju - Jf < 2e-4*Jf, f'certified J exceeds float J by {Ju - Jf:.3e} (relative {(Ju-Jf)/Jf:.2e})')
print(f'time {time.time()-T0:.1f}s')
print(f'{NCHK[0]} checks, {NCHK[1]} failed')
