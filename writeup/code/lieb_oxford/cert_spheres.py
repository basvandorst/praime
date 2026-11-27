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


inp = json.load(open(sys.argv[1]))
RHO_HI = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-3
RHO_LO = float(sys.argv[3]) if len(sys.argv) > 3 else 3e-3
RTAIL = int(sys.argv[4]) if len(sys.argv) > 4 else 200
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
    S2_haa += Haa*d*d/2*i5
    S2_gap += (alpha - (Lm + Lp)/2)*i5
    if k % 1500 == 0:
        print(f'    a = {float(A0):.6f}  G-,G+ = {float(Gm):.8e},{float(Gp):.8e}  Haa+ {float(Haa):.3e}  #sub {len(subs)}  '
              f'iters {i1},{i2}  S2 {float(S2.upper()):.10f}  t {time.time()-t6:.0f}s', flush=True)
print(f'    {len(grid)-1} a-intervals, {nit} bisections, max rel. B&B gap {gapmax:.1e}, {time.time()-t6:.0f}s')
print(f'    contributions: Haa remainder {float(S2_haa.mid()):.3e}, B&B gaps {float(S2_gap.mid()):.3e}')
print('    Haa+ per block:', ' '.join(f'{float(grid[k*BLK]):.3g}:{float(v):.3g}' for k, v in sorted(haa_cache.items())[::3]))
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
print(f'{NCHK[0]} checks, {NCHK[1]} failed')
