#!/usr/bin/env python3
"""Certified upper bounds for the Lieb-Oxford constant (one self-contained verification script).

Part A: mu = nu = uniform unit ball, f = g (LLS 2022, Lemma 10)  ->  c_LO <= 1.604335.
Part B: mu = mixture of 48 spheres, nu = one sphere (exact rational data below), f = g  ->  c_LO <= 1.597856.
Run: python3 lo_certify.py      (about 7 minutes, < 1 GB)
Docstrings of the two parts follow.
"""
MEASURE = {"radii": ["1/24", "1/16", "1/12", "5/48", "1/8", "7/48", "1/6", "3/16", "5/24", "11/48", "1/4", "13/48", "7/24", "5/16", "1/3", "17/48", "3/8", "19/48", "5/12", "7/16", "11/24", "23/48", "1/2", "25/48", "13/24", "9/16", "7/12", "29/48", "5/8", "31/48", "2/3", "11/16", "2839/4096", "17/24", "35/48", "3/4", "37/48", "19/24", "13/16", "5/6", "41/48", "7/8", "43/48", "11/12", "15/16", "23/24", "47/48", "1"], "mu": ["189/1048576", "279/1048576", "267/524288", "921/1048576", "1307/1048576", "1881/1048576", "2435/1048576", "2973/1048576", "3807/1048576", "1149/262144", "5285/1048576", "6501/1048576", "7595/1048576", "8563/1048576", "647/65536", "367/32768", "6535/524288", "1889/131072", "2139/131072", "18825/1048576", "21211/1048576", "2919/131072", "25723/1048576", "16721/524288", "41101/1048576", "26211/524288", "11559/262144", "40063/1048576", "37719/1048576", "2313/65536", "36513/1048576", "36529/1048576", "0", "17173/524288", "2069/65536", "32357/1048576", "32389/1048576", "16007/524288", "3949/131072", "31023/1048576", "15563/524288", "31361/1048576", "31431/1048576", "31381/1048576", "31361/1048576", "30145/1048576", "33259/1048576", "37287/1048576"], "nu": ["0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "1", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0"]}


def part_A():
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
    RHO = 2.5e-4
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
    DUMP_A = []
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
        DUMP_A.append((A0, A1, float(Gm.upper()), float(Gp.upper()), float(Haa.upper())))
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
    import json as _js; _js.dump(DUMP_A, open('dump_A.json', 'w'))
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
    print(f'part A: {NCHK[0]} checks, {NCHK[1]} failed')
    return NCHK


def part_B():
    #!/usr/bin/env python3
    """Certified upper bound for the Lieb-Oxford constant with mu, nu finite mixtures of uniform sphere measures
    (exact rational radii and weights, read from a JSON file) and f = g of LLS 2022, Lemma 10, eq. (51).

       c_LO <= (3/2) (2 D(mu,mu)^2 J)^(1/3),  J = 4 pi int_0^oo g(r) r^-5 dr.

    Building block: the sphere formula (43), proved exactly below (sympy).  For spheres of radius r (at 0, scale a)
    and q (at e1, scale b), T = r^3 q^3 Phi_ss(a/r, b/q) is the polynomial
       I   a^3 b^3 (1 - b/q)                 if  aq - br - ab >= 0   (sphere at 0 inside the one at e1)
       II  a^3 b^3 (1 - a/r)                 if  br - aq - ab >= 0   (sphere at e1 inside the one at 0)
       III a^2 b^2 (aq + br - ab)^2/(4 r q)  if  both <= 0 <= X1 = aq + br - ab
       0                                     if  X1 <= 0             (disjoint)
    Psi_{mu nu} = sum over atom pairs c_{rq} T_{rq},  c_{rq} = mu_r nu_q + nu_r mu_q - nu_r nu_q  (eq. (45), (69)).
    Run: python3 cert_spheres.py measure.json [rho_hi rho_lo R]
    """

    import sys, json, time, math, heapq, bisect
    from fractions import Fraction
    import sympy as sp
    from flint import arb, arb_poly, fmpq, ctx

    ctx.prec = 128
    T0 = time.time()
    NCHK = [0, 0]


    def check(cond, msg):
        NCHK[0] += 1
        print(('ok   ' if cond else 'FAIL ') + msg, flush=True)
        if not cond:
            NCHK[1] += 1


    def Q(x):
        x = Fraction(x)
        return fmpq(x.numerator, x.denominator)


    inp = MEASURE
    RHO_HI = 1e-3
    RHO_LO = 3e-3
    RTAIL = 600
    radF = [Fraction(x) for x in inp['radii']]
    muF = [Fraction(x) for x in inp['mu']]
    nuF = [Fraction(x) for x in inp['nu']]
    print('atoms', [str(r) for r in radF])
    check(sum(muF) == 1 and sum(nuF) == 1 and min(muF) >= 0 and min(nuF) >= 0, 'mu, nu are probability vectors (exact)')
    check(all(0 < r <= 1 for r in radF) and len(set(radF)) == len(radF), 'radii distinct, in (0,1]: supports in the unit ball, mu({0}) = 0')

    # ---------------- Part 1: exact proof of (43) and of the piece formulas ----------------
    print('--- Part 1: sphere formula (43) and the pieces')
    a, b, s, r_, q_ = sp.symbols('a b s r q', positive=True)
    R1, R2 = sp.symbols('R1 R2', positive=True)
    ps = s/(2*R2)                         # density of |y| for y uniform on the sphere of radius R2 about e1
    PHI = lambda W: (1 - W)/(R1**3*R2**3)
    A_, B_ = 1/R1, 1/R2
    sig = lambda amb: A_**2*B_**2/4*((A_ + B_ - A_*B_)**2 - (amb - A_*B_)**2)
    one = A_**2*B_**2/4*(A_ + B_ - A_*B_)**2
    tests = [('disjoint', PHI(sp.integrate(ps/s, (s, 1 - R2, 1 + R2))), 0),
             ('1 inside 2', PHI(sp.integrate(ps/s, (s, R2 - 1, 1 + R2))), sig(A_ - B_)),
             ('2 inside 1, R2<=1', PHI(sp.integrate(ps/R1, (s, 1 - R2, 1 + R2))), sig(B_ - A_)),
             ('2 inside 1, R2>=1', PHI(sp.integrate(ps/R1, (s, R2 - 1, 1 + R2))), sig(B_ - A_)),
             ('overlap, R2<=1', PHI(sp.integrate(ps/R1, (s, 1 - R2, R1)) + 0*s + sp.integrate(ps/s, (s, R1, 1 + R2))) - 0, one),
             ('overlap, R2>=1', PHI(sp.integrate(ps/R1, (s, R2 - 1, R1)) + sp.integrate(ps/s, (s, R1, 1 + R2))), one)]
    for name, lhs, rhs in tests:
        check(sp.simplify(sp.expand(lhs - rhs)) == 0, f'(43) in configuration {name}')
    x, y = a/r_, b/q_
    X1, X2p, X2m = a*q_ + b*r_ - a*b, a*q_ - b*r_ - a*b, b*r_ - a*q_ - a*b
    base = r_**3*q_**3*x**2*y**2/4
    TI, TII = a**3*b**3*(1 - b/q_), a**3*b**3*(1 - a/r_)
    TIII = a**2*b**2*X1**2/(4*r_*q_)
    check(sp.simplify(base*((X1/(r_*q_))**2 - (X2p/(r_*q_))**2) - TI) == 0, 'piece I = r^3q^3 Phi_ss when x2 = (aq-br-ab)/(rq) >= 0')
    check(sp.simplify(base*((X1/(r_*q_))**2 - (X2m/(r_*q_))**2) - TII) == 0, 'piece II = r^3q^3 Phi_ss when (br-aq-ab)/(rq) >= 0')
    check(sp.simplify(base*(X1/(r_*q_))**2 - TIII) == 0, 'piece III = r^3q^3 Phi_ss in the overlap')
    # hand-coded coefficients in b (degree <= 4) of T, T_a, T_aa for each piece, checked against sympy
    def coef_expr(piece, order):
        e = {'I': TI, 'II': TII, 'III': TIII}[piece]
        return sp.expand(sp.diff(e, a, order))
    HAND = {
        ('I', 0): lambda A, r, q: [0, 0, 0, A**3, -A**3/q],
        ('I', 1): lambda A, r, q: [0, 0, 0, 3*A**2, -3*A**2/q],
        ('I', 2): lambda A, r, q: [0, 0, 0, 6*A, -6*A/q],
        ('II', 0): lambda A, r, q: [0, 0, 0, A**3 - A**4/r, 0],
        ('II', 1): lambda A, r, q: [0, 0, 0, 3*A**2 - 4*A**3/r, 0],
        ('II', 2): lambda A, r, q: [0, 0, 0, 6*A - 12*A**2/r, 0],
        ('III', 0): lambda A, r, q: [0, 0, A**4*q**2/(4*r*q), 2*A**3*q*(r - A)/(4*r*q), A**2*(r - A)**2/(4*r*q)],
        ('III', 1): lambda A, r, q: [0, 0, 4*A**3*q**2/(4*r*q), (6*A**2*q*r - 8*A**3*q)/(4*r*q), (2*A*r**2 - 6*A**2*r + 4*A**3)/(4*r*q)],
        ('III', 2): lambda A, r, q: [0, 0, 12*A**2*q**2/(4*r*q), (12*A*q*r - 24*A**2*q)/(4*r*q), (2*r**2 - 12*A*r + 12*A**2)/(4*r*q)],
    }
    okh = True
    for (pc, od), fn in HAND.items():
        ex = coef_expr(pc, od)
        hand = sum(c*b**j for j, c in enumerate(fn(a, r_, q_)))
        okh = okh and sp.simplify(sp.expand(ex - hand)) == 0
    check(okh, 'hand-coded b-coefficients of T, T_a, T_aa agree with sympy (all pieces)')
    uu = q_ - b
    check(sp.expand(b**2*(12*a**2*uu**2 + 12*a*uu*b*r_ + 2*b**2*r_**2)/(4*r_*q_) - sp.diff(TIII, a, 2)) == 0,
          'factored form of d^2 T_III/da^2')

    # ---------------- Part 2: measure data ----------------
    n = len(radF)
    pairs = []   # (r, q, c) exact
    for i in range(n):
        for k in range(n):
            c = muF[i]*nuF[k] + nuF[i]*muF[k] - nuF[i]*nuF[k]
            if c != 0:
                pairs.append((radF[i], radF[k], c))
    print(f'    {n} atoms, {len(pairs)} pairs with c != 0')
    check(sum(c for _, _, c in pairs) == 1, 'sum of pair coefficients = 1 (Psi ~ a^3 b^3 at the origin)')
    DmuF = sum(muF[i]*muF[k]/(2*max(radF[i], radF[k])) for i in range(n) for k in range(n))
    print(f'    D(mu,mu) = {DmuF} = {float(DmuF):.12f}')
    pairsA = [(Q(r), Q(qq), Q(c), arb(Q(r)), arb(Q(qq)), arb(Q(c))) for (r, qq, c) in pairs]


    def coefs(pc, od, A, ra, qa):
        return HAND[(pc, od)](A, ra, qa)


    # diagonal Psi(b,b) = sum c T(b,b):  I b^6(1-b/q) if b <= q-r; II b^6(1-b/r) if b <= r-q;
    # III b^6 (q+r-b)^2/(4rq) if |q-r| <= b <= q+r; 0 beyond.
    diag_bp = sorted(set([Fraction(0)] + [abs(qq - r) for r, qq, _ in pairs if qq != r] + [qq + r for r, qq, _ in pairs]))
    bsym = sp.symbols('bb')
    diag_polys = []   # list of (lo, hi, arb_poly) on [lo, hi]; last segment [max, oo) is 0
    for lo, hi in zip(diag_bp[:-1], diag_bp[1:]):
        m = (lo + hi)/2
        e = 0
        for r, qq, c in pairs:
            if m <= qq - r:
                e += c*bsym**6*(1 - bsym/qq)
            elif m <= r - qq:
                e += c*bsym**6*(1 - bsym/r)
            elif m <= qq + r:
                e += c*bsym**6*(qq + r - bsym)**2/(4*r*qq)
        P = sp.Poly(sp.expand(e), bsym)
        diag_polys.append((lo, hi, arb_poly([arb(Q(sp.Rational(cc))) for cc in reversed(P.all_coeffs())])))
    DIAG_MAX = diag_bp[-1]
    diag_lo = [d[0] for d in diag_polys]
    ZERO8 = arb_poly([0])


    def diag_poly_at(mid):
        """polynomial of Psi(b,b) on the diagonal segment containing the (exact) point mid"""
        if mid >= DIAG_MAX:
            return ZERO8, None
        i = bisect.bisect_right(diag_lo, mid) - 1
        return diag_polys[i][2], (diag_polys[i][0], diag_polys[i][1])


    # ---------------- Part 3: branch and bound in b for fixed a = am ----------------
    def subintervals(amF, bmax):
        """exact breakpoints in b at a = am (Fractions) and, for each sub-interval, the Psi and Psi_a b-polynomials"""
        am = Q(amF)
        A = arb(am)
        events = []
        for idx, (r, qq, c) in enumerate(pairs):
            events.append((amF*qq/(amF + r), idx, 'III'))
            if amF < r:
                events.append((amF*qq/(r - amF), idx, 'II'))
            elif amF > r:
                events.append((amF*qq/(amF - r), idx, '0'))
        events.sort(key=lambda e: e[0])
        cuts = sorted(set([Fraction(0), bmax] + [e[0] for e in events if 0 < e[0] < bmax] +
                          [d for d in diag_bp if 0 < d < bmax]))
        # per-pair coefficient lists at A
        cf = []
        for (r, qq, c, ra, qa, ca) in pairsA:
            cf.append({pc: ([ca*x for x in coefs(pc, 0, A, ra, qa)], [ca*x for x in coefs(pc, 1, A, ra, qa)])
                       for pc in ('I', 'II', 'III')})
        state = ['I']*len(pairs)
        S0 = [arb(0)]*5
        S1 = [arb(0)]*5
        for idx in range(len(pairs)):
            S0 = [u + v for u, v in zip(S0, cf[idx]['I'][0])]
            S1 = [u + v for u, v in zip(S1, cf[idx]['I'][1])]
        out = []
        ev = 0
        for lo, hi in zip(cuts[:-1], cuts[1:]):
            # apply all events with breakpoint <= lo (pieces are constant on (lo, hi))
            while ev < len(events) and events[ev][0] <= lo:
                _, idx, new = events[ev]
                old = state[idx]
                if old != '0':
                    S0 = [u - v for u, v in zip(S0, cf[idx][old][0])]
                    S1 = [u - v for u, v in zip(S1, cf[idx][old][1])]
                if new != '0':
                    S0 = [u + v for u, v in zip(S0, cf[idx][new][0])]
                    S1 = [u + v for u, v in zip(S1, cf[idx][new][1])]
                state[idx] = new
                ev += 1
            out.append((lo, hi, arb_poly(S0), arb_poly(S1)))
        return out


    def piece_check(amF, lo, hi, P0):
        """independent recomputation of Psi(am, mid) by direct piece selection (used for spot checks)"""
        m = (lo + hi)/2
        v = Fraction(0)
        for r, qq, c in pairs:
            x1, x2p, x2m = amF*qq + m*r - amF*m, amF*qq - m*r - amF*m, m*r - amF*qq - amF*m
            if x2p >= 0:
                v += c*amF**3*m**3*(1 - m/qq)
            elif x2m >= 0:
                v += c*amF**3*m**3*(1 - amF/r)
            elif x1 >= 0:
                v += c*amF**2*m**2*x1**2/(4*r*qq)
        return abs(float(P0(arb(Q(m))).mid()) - float(v)) <= 1e-12*max(1.0, abs(float(v)))


    def taylor_ub(Pc, P1c, P2box, rad):
        u2 = P2box.upper()
        return (Pc.upper() + abs(P1c).upper()*rad + (u2 if u2 > 0 else arb(0))*(rad*rad/2)).upper()


    def bb_max(subs, sv, tol_rel=1e-10, maxit=40000):
        """certified upper bound of max_b F(b), F = Psi(am,b) + sv Psi_a(am,b) - Psi(b,b)_+/2, b in [0,bmax]"""
        heap = []
        L = [None]
        polys = []
        for (lo, hi, P0, P1) in subs:
            P = P0 + P1*sv
            Pd, seg = diag_poly_at((lo + hi)/2)
            polys.append((P, P.derivative(), P.derivative().derivative(), Pd, Pd.derivative(), Pd.derivative().derivative()))

        def push(i, lo, hi):
            P, P1_, P2_, Pd, Pd1, Pd2 = polys[i]
            c = arb(Q((lo + hi)/2))
            rad = arb(Q((hi - lo)/2)).upper()
            box = arb(Q(lo)).union(arb(Q(hi)))
            Pc, P1c, P2b = P(c), P1_(c), P2_(box)
            Dc, D1c, D2b = Pd(c), Pd1(c), Pd2(box)
            ubA = taylor_ub(Pc, P1c, P2b, rad)                                     # F <= P
            ubB = taylor_ub(Pc - Dc/2, P1c - D1c/2, P2b - D2b/2, rad)             # F <= P - Pd/2
            U = ubA if ubA < ubB else ubB
            lowc = Pc.lower() - (Dc.upper() if Dc.upper() > 0 else arb(0))/2
            if L[0] is None or lowc > L[0]:
                L[0] = lowc
            heapq.heappush(heap, (-float(U), float(lo), lo, hi, i, U))
        for i, (lo, hi, _, _) in enumerate(subs):
            push(i, lo, hi)
        it = 0
        while True:
            top = heap[0]
            U, lo, hi, i = top[5], top[2], top[3], top[4]
            tol = tol_rel*max(abs(float(L[0])), 1e-300)
            if float(U) <= float(L[0]) + tol or it >= maxit or (hi - lo) < Fraction(1, 2**60):
                best = U
                for e in heap:
                    if e[5] > best:
                        best = e[5]
                return best, it, float(L[0])
            heapq.heappop(heap)
            m = (lo + hi)/2
            # bisect at a short dyadic-ish point to keep fractions small
            m = Fraction(round(float(m)*2**52), 2**52) if lo < Fraction(round(float(m)*2**52), 2**52) < hi else m
            push(i, lo, m)
            push(i, m, hi)
            it += 1


    # ---------------- Part 4: Haa+ on blocks ----------------
    def taa_box(A, B, corners, r, qq, ra, qa):
        """enclosure of T_aa on the box A x B.  Piece tests: X1, X2p, X2m are bilinear in (a,b), so their extreme values
        on the box are at the 4 corners; a piece is dropped only if excluded at all corners with a safety margin
        (keeping a piece is always safe: the hull only grows)."""
        rf, qf = float(r), float(qq)
        v1 = [x*qf + y*rf - x*y for x, y in corners]
        v2p = [x*qf - y*rf - x*y for x, y in corners]
        v2m = [y*rf - x*qf - x*y for x, y in corners]
        eps = 1e-9*max((1 + x)*(1 + y) for x, y in corners)
        vals = []
        if max(v2p) > -eps:
            vals.append(6*A*B**3 - 6*A*B**4/qa)
        if max(v2m) > -eps:
            vals.append((6*A - 12*A*A/ra)*B**3)
        if max(v1) > -eps and min(v2p) < eps and min(v2m) < eps:
            U = qa - B          # factored form: T_III,aa = b^2 (12 a^2 u^2 + 12 a u b r + 2 b^2 r^2)/(4 r q), u = q - b
            vals.append(B*B*(12*A*A*U*U + 12*A*U*B*ra + 2*B*B*ra*ra)/(4*ra*qa))
        if min(v1) < eps:
            vals.append(arb(0))
        v = vals[0]
        for w in vals[1:]:
            v = v.union(w)
        return v


    RMAX = max(radF)


    def haa_sup(P, Qe, nb=160, na=1):
        """upper bound of max(0, ess sup Psi_aa) on [P,Qe] x [0,Qe]; P, Qe Fractions"""
        bsup = Qe
        if P > RMAX:
            bsup = min(Qe, max(P*qq/(P - r) for r, qq, _ in pairs))   # Psi = 0 beyond (x1 <= 0 for every pair)
        best = arb(0)
        for i in range(na):
            alo, ahi = P + (Qe - P)*i/na, P + (Qe - P)*(i + 1)/na
            A = arb(Q(alo)).union(arb(Q(ahi)))
            for j in range(nb):
                blo, bhi = bsup*j/nb, bsup*(j + 1)/nb
                B = arb(Q(blo)).union(arb(Q(bhi)))
                corners = [(float(x), float(y)) for x in (alo, ahi) for y in (blo, bhi)]
                tot = arb(0)
                for (r, qq, c, ra, qa, ca) in pairsA:
                    tot += ca*taa_box(A, B, corners, r, qq, ra, qa)
                u = tot.upper()
                if u > best:
                    best = u
        return best


    # ---------------- Part 5: S2 ----------------
    print('--- S2: Taylor in a, convex chord, branch and bound in b')
    A0F = Fraction(1, 4096)
    grid = [A0F]
    while grid[-1] < RTAIL:
        x_ = grid[-1]
        rho = RHO_LO if x_ < Fraction(1, 4) else RHO_HI
        nxt = Fraction(float(x_)*(1 + rho))
        grid.append(min(nxt, Fraction(RTAIL)))
    check(grid[0] == A0F and grid[-1] == RTAIL and all(u < v for u, v in zip(grid[:-1], grid[1:])), f'a-grid covers [1/4096, {RTAIL}] ({len(grid)-1} intervals)')
    BLK = 32
    DUMP_B = []
    haa_cache = {}
    S2 = arb(0)
    t6 = time.time()
    nit = 0
    gapmax = 0.0
    spot_ok = True
    S2_haa = arb(0)
    S2_gap = arb(0)
    for k in range(len(grid) - 1):
        A0, A1 = grid[k], grid[k + 1]
        blk = k//BLK
        if blk not in haa_cache:
            haa_cache[blk] = haa_sup(grid[blk*BLK], grid[min((blk + 1)*BLK, len(grid) - 1)])
        Haa = haa_cache[blk]
        amF = (A0 + A1)/2
        dF = A1 - amF
        subs = subintervals(amF, A1)
        if k % 997 == 0:
            for (lo, hi, P0, P1) in subs[::max(1, len(subs)//7)]:
                spot_ok = spot_ok and piece_check(amF, lo, hi, P0)
        d = arb(Q(dF))
        Gm, i1, Lm = bb_max(subs, -d)
        Gp, i2, Lp = bb_max(subs, d)
        nit += i1 + i2
        gapmax = max(gapmax, (float(Gm) - Lm)/max(abs(Lm), 1e-300), (float(Gp) - Lp)/max(abs(Lp), 1e-300))
        alpha, beta = (Gm + Gp)/2, (Gp - Gm)/(2*d)
        A0b, A1b, amb = arb(Q(A0)), arb(Q(A1)), arb(Q(amF))
        i5 = (A0b**-4 - A1b**-4)/4
        i4 = (A0b**-3 - A1b**-3)/3
        S2 += alpha*i5 + beta*(i4 - amb*i5) + Haa*d*d/2*i5
        DUMP_B.append((float(A0), float(A1), str(A0), str(A1), float(Gm.upper()), float(Gp.upper()), float(Haa.upper()), blk))
        S2_haa += Haa*d*d/2*i5
        S2_gap += (alpha - (Lm + Lp)/2)*i5
        if k % 1500 == 0:
            print(f'    a = {float(A0):.6f}  G-,G+ = {float(Gm):.8e},{float(Gp):.8e}  Haa+ {float(Haa):.3e}  #sub {len(subs)}  '
                  f'iters {i1},{i2}  S2 {float(S2.upper()):.10f}  t {time.time()-t6:.0f}s', flush=True)
    print(f'    {len(grid)-1} a-intervals, {nit} bisections, max rel. B&B gap {gapmax:.1e}, {time.time()-t6:.0f}s')
    print(f'    contributions: Haa remainder {float(S2_haa.mid()):.3e}, B&B gaps {float(S2_gap.mid()):.3e}')
    print('    Haa+ per block:', ' '.join(f'{float(grid[k*BLK]):.3g}:{float(v):.3g}' for k, v in sorted(haa_cache.items())[::3]))
    import json as _js; _js.dump(DUMP_B, open('dump_B.json', 'w'))
    check(spot_ok, 'spot checks: event-built Psi polynomials = direct exact piece evaluation')
    # sanity checks by random float sampling (not part of the proof)
    import random
    random.seed(12345)
    pairs_f = [(float(r), float(qq), float(c)) for r, qq, c in pairs]


    def term_f(x, y, r, qq, order):
        """T (order 0), T_a (1), T_aa (2) of one pair, float, piece selected at (x,y)"""
        x1, x2p, x2m = x*qq + y*r - x*y, x*qq - y*r - x*y, y*r - x*qq - x*y
        if x2p >= 0:
            pc = 'I'
        elif x2m >= 0:
            pc = 'II'
        elif x1 >= 0:
            pc = 'III'
        else:
            return 0.0
        return sum(cf*y**j for j, cf in enumerate(HAND[(pc, order)](x, r, qq)))


    def psi_f(x, y, order):
        return sum(c*term_f(x, y, r, qq, order) for r, qq, c in pairs_f)


    def diag_f(y):
        return sum(c*term_f(y, y, r, qq, 0) for r, qq, c in pairs_f)


    bad = 0
    for blk, Hv in haa_cache.items():
        P = float(grid[blk*BLK]); Qv = float(grid[min((blk + 1)*BLK, len(grid) - 1)])
        for _ in range(10):
            x, y = random.uniform(P, Qv), random.uniform(0, min(Qv, 2.5))
            if psi_f(x, y, 2) > float(Hv)*(1 + 1e-9) + 1e-12:
                bad += 1
    check(bad == 0, f'sampling: Psi_aa <= Haa+ at {10*len(haa_cache)} random points of the blocks')
    bad = 0
    for k in range(0, len(grid) - 1, 1499):
        A0, A1 = grid[k], grid[k + 1]
        amF = (A0 + A1)/2; dF = A1 - amF
        subs = subintervals(amF, A1)
        for sgn in (-1, 1):
            Gv, _, _ = bb_max(subs, arb(Q(dF))*sgn)
            am, dd = float(amF), float(dF)
            top = min(float(A1), 2.5)
            fm = max(psi_f(am, y, 0) + sgn*dd*psi_f(am, y, 1) - max(diag_f(y), 0)/2 for y in [top*j/1500 for j in range(1501)])
            if fm > float(Gv) + 1e-12*abs(fm):
                bad += 1
    check(bad == 0, 'sampling: certified G(+-d) >= float maximum on a 1501-point b-grid (every 1499th interval)')

    # ---------------- Part 6: S0 and tail ----------------
    print('--- S0, tail')
    # a <= a0: Psi <= Phi_mu_nu + Phi_nu_mu <= 2 a^3 b^3 (Phi_nu_nu >= 0 by (21); Phi <= a^3 b^3), so g(a) <= 2 a^6
    S0 = arb(Q(A0F))**2           # int_0^a0 2 r dr
    # kappa_mu = sup_b Lambda(b), Lambda(b) = b^3 (1 - sum_k mu_k min(1, b/rho_k)); Lambda = 0 for b >= max rho
    mu_atoms = sorted((r, m) for r, m in zip(radF, muF) if m > 0)
    segs = sorted(set([Fraction(0)] + [r for r, _ in mu_atoms]))
    lam_polys = []
    for lo, hi in zip(segs[:-1], segs[1:]):
        m = (lo + hi)/2
        inner = sum(w for r, w in mu_atoms if r <= m)
        slope = sum(w/r for r, w in mu_atoms if r > m)
        # Lambda = b^3 (1 - inner - slope b) = (1-inner) b^3 - slope b^4
        P = arb_poly([0, 0, 0, arb(Q(1 - inner)), -arb(Q(slope))])
        lam_polys.append((lo, hi, P))


    def lam_max(lo_w=Fraction(0), hi_w=None):
        """certified upper bound of sup Lambda on [lo_w, hi_w] (whole support if hi_w is None)"""
        heap = []
        L = [None]

        def push(P, lo, hi):
            c = arb(Q((lo + hi)/2)); rad = arb(Q((hi - lo)/2)).upper(); box = arb(Q(lo)).union(arb(Q(hi)))
            U = taylor_ub(P(c), P.derivative()(c), P.derivative().derivative()(box), rad)
            lw = P(c).lower()
            if L[0] is None or lw > L[0]:
                L[0] = lw
            heapq.heappush(heap, (-float(U), float(lo), lo, hi, P, U))
        for lo, hi, P in lam_polys:
            lo2, hi2 = max(lo, lo_w), (hi if hi_w is None else min(hi, hi_w))
            if lo2 < hi2:
                push(P, lo2, hi2)
        if not heap:
            return arb(0), 0.0
        while True:
            top = heap[0]
            if float(top[5]) <= float(L[0]) + 1e-13*abs(float(L[0])) or top[3] - top[2] < Fraction(1, 2**60):
                best = top[5]
                for e in heap:
                    if e[5] > best:
                        best = e[5]
                return best, float(L[0])
            heapq.heappop(heap)
            m = (top[2] + top[3])/2
            push(top[4], top[2], m)
            push(top[4], m, top[3])


    kap, kapL = lam_max()
    print(f'    kappa_mu in [{kapL:.15f}, {float(kap):.15f}]')
    # windows: for a >= R and pair (r,q), Psi-term/a^3 differs from c L_q(b) only for b in [q/(1+r/R), q/(1-r/R)],
    # by at most |c| r (q/(1-r/R))^3/(1-r/R) / a.  Group by q, merge overlapping windows.
    check(RTAIL > 2*RMAX, 'R > 2 max radius (region II impossible, windows well defined)')
    tR = Fraction(1, RTAIL)
    groups = {}
    for r, qq, c in pairs:
        groups.setdefault(qq, []).append((r, c))
    wins = []
    for qq, lst in groups.items():
        rmx = max(r for r, _ in lst)
        wins.append([qq/(1 + rmx*tR), qq/(1 - rmx*tR), sum(abs(c)*r*(qq/(1 - r*tR))**3/(1 - r*tR) for r, c in lst)])
    wins.sort()
    clusters = []
    for w in wins:
        if clusters and w[0] <= clusters[-1][1]:
            clusters[-1][1] = max(clusters[-1][1], w[1]); clusters[-1][2] += w[2]
        else:
            clusters.append(list(w))
    MU_MAXR = max(r for r, m in mu_atoms)
    M_R = kap
    worst = None
    for lo, hi, C in clusters:
        sw, _ = lam_max(lo, hi)
        if hi > MU_MAXR and sw < 0:
            sw = arb(0)          # Lambda = 0 beyond the support of mu
        U = (sw + arb(Q(C))/RTAIL).upper()
        if U > M_R:
            M_R, worst = U, (float(lo), float(hi), float(C))
    S3 = M_R/arb(RTAIL)
    print(f'    {len(clusters)} window clusters; M_R = {float(M_R):.15f} (kappa {float(kap):.15f}; worst cluster {worst}); tail <= M_R/R = {float(S3.upper()):.12e}')

    # ---------------- Part 7: assembly ----------------
    J = 4*arb.pi()*(S0 + S2 + S3)
    D = arb(Q(DmuF))
    cLO = arb(Q(Fraction(3, 2)))*(2*J*D*D).root(3)
    print(f'    S0 <= {float(S0.upper()):.3e}, S2 <= {float(S2.upper()):.12f}, S3 <= {float(S3.upper()):.12f}')
    print(f'    J <= {J.upper()}')
    print(f'    D(mu,mu) = {float(DmuF):.15f}')
    print(f'    c_LO <= {cLO.upper()}')
    cu = math.ceil(float(cLO.upper())*1e6)/1e6
    check(cu >= float(cLO.upper()), f'outward rounding: c_LO <= {cu:.6f}')

    # float cross-check (independent numpy code, dense grids)
    import numpy as np
    rr = np.array([float(x) for x in radF]); mm = np.array([float(x) for x in muF]); nn = np.array([float(x) for x in nuF])
    cc = np.outer(mm, nn) + np.outer(nn, mm) - np.outer(nn, nn)


    def psi_np(A, B):
        out = np.zeros(np.broadcast(A, B).shape)
        for i in range(len(rr)):
            for k in range(len(rr)):
                if cc[i, k] != 0:
                    X1 = A*rr[k] + B*rr[i] - A*B
                    X2 = np.abs(A*rr[k] - B*rr[i]) - A*B
                    out += cc[i, k]*A*A*B*B/(4*rr[i]*rr[k])*(np.maximum(X1, 0)**2 - np.maximum(X2, 0)**2)
        return out


    xa = np.concatenate([np.linspace(1e-4, 4, 4000), np.exp(np.linspace(np.log(4), np.log(RTAIL), 800))[1:]])
    gv = np.empty_like(xa)
    for i, x_ in enumerate(xa):
        bsupf = x_ if x_ <= 1.01*rr.max() else min(x_, max(x_*rr[k]/(x_ - rr[i]) for i in range(len(rr)) for k in range(len(rr)) if cc[i, k] != 0))
        bb_ = np.linspace(0, bsupf, 3001)
        h = psi_np(np.full_like(bb_, x_), bb_) - np.maximum(psi_np(bb_, bb_), 0)/2
        gv[i] = h.max()
    Jf = 4*np.pi*(np.trapezoid(gv/xa**5, xa) + gv[0]/xa[0]**5*xa[0]/2 + kapL/RTAIL)
    cf = 1.5*(2*Jf*float(DmuF)**2)**(1/3)
    print(f'    float cross-check (grid): J = {Jf:.8f}, c = {cf:.8f}')
    check(Jf <= float(J.upper()) and float(J.upper()) - Jf < 1e-3*Jf, f'certified J above float J by {(float(J.upper())-Jf)/Jf:.2e} (relative)')
    print(f'time {time.time()-T0:.0f}s')
    print(f'part B: {NCHK[0]} checks, {NCHK[1]} failed')
    return NCHK


if __name__ == '__main__':
    nA = part_A()
    nB = part_B()
    print(f'{nA[0] + nB[0]} checks, {nA[1] + nB[1]} failed')
