"""Rigorous certificate (python-flint / arb ball arithmetic) for one class-A solution.

Input: polish_A_d<d>.json (KKT optimum, b_j to 40 digits). We round b_j to 25 significant digits (exact rationals)
and add eps to b_0 (g -> g + eps*Phi*Phi, B -> B + eps). Then we PROVE (ball arithmetic, own Bessel K: kbessel.py):
  (1) B(t) = sum b_j t^{2j} > 0 for all real t      [exact: Sturm count (sympy) of roots of q(s)=sum b_j s^j in [0,oo)]
  (2) g(u) > 0 on [0, U]       [subdivision; on [c-r, c+r]: Taylor polynomial of order p at c with ball derivatives
                                D^k g(c), k<=p; exact minimum of its quadratic part; |terms k=3..p| <= |D^k g(c)| r^k/k!;
                                remainder <= sup|D^{p+1} g| r^{p+1}/(p+1)!  (crude: monotonicity of K0, K1 and of |P|)]
  (3) g(u) > 0 on [U, oo)      [asymptotic expansion of K0, K1 with the DLMF 10.40.10 remainder, N >= 2 terms bounded]
  (4) cost(g) = 2 sum_n Lambda(n) n^{-1/2} g(log n) <= C_up,  g(0) >= g0_lo       [balls + explicit tail bounds]
Representation used (exact, see README): g(u) = sum_N T_N(u),
  T_N(u) = 4 d(N) (pi N)^{-1/2} X^{1/2} [P0(X) K0(2X) + P1(X) K1(2X)],  X = pi N e^u,
P0, P1 rational polynomials of g = sum_j b_j (-1)^j D^{2j}(Phi*Phi), D = d/du = X d/dX  (core.basis_polys).
Bound for |T_N| (u >= 0): K0 <= K1 <= sqrt(pi/(2x)) e^{-x} (1 + 3/(8x)) (DLMF 10.40.10, l=1), d(N) <= 2 sqrt(N):
  |T_N| <= h(N) = 4.4 Pabs(1) X^M e^{-2X},  Pabs = |P0|+|P1| coefficientwise, M = max degree (X >= 1); for N >= M
  the ratios h(N+1)/h(N) are <= e^{1-2pi} < 0.0051; tails are summed explicitly up to N = M, then geometrically.
Usage: python certify.py d [eps] [U] [p]"""
import sys, json, time, os
from fractions import Fraction
from flint import arb, arb_poly, fmpq, ctx
import sympy
from core import basis_polys, divisor_count, von_mangoldt
from kbessel import K01      # own rigorous K0, K1 (arb's bessel_k gave a useless ball at x ~ 153, see arbrepro.py)

ctx.prec = int(os.environ.get("PREC", 800))   # >= target + 2.9*120 bits: the series branch of K01 amplifies the input radius by e^{2x}
d = int(sys.argv[1])
EPS = Fraction(sys.argv[2]) if len(sys.argv) > 2 else Fraction(1, 10 ** 6)
U = Fraction(sys.argv[3]) if len(sys.argv) > 3 else Fraction(3)
PT = int(sys.argv[4]) if len(sys.argv) > 4 else 12           # Taylor order
SRC = "polish_A_d%d.json" % d if os.path.exists("polish_A_d%d.json" % d) else "results_A_d%d.json" % d
R = json.load(open(SRC))
SIG = int(os.environ.get("SIG", 25))     # significant digits kept from the 40-digit b_j
b = [Fraction(str(sympy.Float(sympy.Float(s, max(45, SIG + 5)), SIG))) for s in R["b"]]
b[0] += EPS
n = d + 1
log = dict(d=d, deg=2 * d, source=SRC, eps=str(EPS), U=str(U), taylor_order=PT, significant_digits=SIG, b_rational=[str(x) for x in b])
t0 = time.time()


def fq(x):
    return arb(fmpq(x.numerator, x.denominator))


# ---------------- (1) B > 0 on R, exact ----------------
s = sympy.symbols("s")
q = sympy.Poly([sympy.Rational(x.numerator, x.denominator) for x in reversed(b)], s, domain="QQ")
nroots = q.count_roots(0, None)
log["B_roots_in_s>=0"] = int(nroots)
assert nroots == 0 and b[0] > 0, "B not positive"
print("(1) q(s) = sum b_j s^j has %d roots in [0,oo) and q(0) = %.6g > 0, so B(t) = q(t^2) > 0 on R  [exact, Sturm]"
      % (nroots, float(b[0])), flush=True)

# ---------------- polynomials of D^r g, r = 0..PT+1 ----------------
RMAX = PT + 1
BP = basis_polys(d, RMAX)


def strip(c):
    c = list(c)
    while len(c) > 1 and c[-1] == 0:
        c.pop()
    return c


P = []
for r in range(RMAX + 1):
    L0 = max(len(BP[j][r][0]) for j in range(n))
    L1 = max(len(BP[j][r][1]) for j in range(n))
    P0 = [sum((b[j] * BP[j][r][0][k] for j in range(n) if k < len(BP[j][r][0])), Fraction(0)) for k in range(L0)]
    P1 = [sum((b[j] * BP[j][r][1][k] for j in range(n) if k < len(BP[j][r][1])), Fraction(0)) for k in range(L1)]
    P.append((strip(P0), strip(P1)))


def apoly(c):
    return arb_poly([fq(x) for x in c])


AP = [(apoly(P0), apoly(P1)) for (P0, P1) in P]
APabs = [(apoly([abs(x) for x in P0]), apoly([abs(x) for x in P1])) for (P0, P1) in P]
Mdeg = [max(len(P0), len(P1)) - 1 for (P0, P1) in P]
Pabs1 = [fq(sum(abs(x) for x in P0) + sum(abs(x) for x in P1)) for (P0, P1) in P]
PI = arb.pi()
RATIO = arb("0.0051")          # >= e^{1-2pi} = 0.00509...
assert (RATIO - (1 - 2 * PI).exp()).lower() > 0
FACT = [arb(1)]
for k in range(1, RMAX + 2):
    FACT.append(FACT[-1] * k)


def tail_bound(r, a, bnd, N0):
    """bound for sum_{N>N0} |T_N^{(r)}(u)|, valid for all u in [a, bnd] (a >= 0):
    sum_{N=N0+1}^{N1-1} h(N) + h(N1)/(1-0.0051),  N1 = max(N0, M)+1,  h(N) = 4.4 Pabs(1) (pi N e^bnd)^M e^{-2 pi N e^a}."""
    M = Mdeg[r]
    ea = fq(a).exp()
    eb = fq(bnd).exp()
    N1 = max(N0, M) + 1
    tot = arb(0)
    for N in range(N0 + 1, N1 + 1):
        h = arb("4.4") * Pabs1[r] * (PI * N * eb) ** M * (-2 * PI * N * ea).exp()
        tot += h / (1 - RATIO) if N == N1 else h
    return tot.upper()


def point_vals(u, rs):
    """balls for D^r g(u), r in rs (u exact rational >= 0)."""
    eu = fq(u).exp()
    out = {r: arb(0) for r in rs}
    M = max(Mdeg[r] for r in rs)
    N = 1
    while True:
        X = PI * N * eu
        k0, k1 = K01(2 * X)
        pref = 4 * divisor_count(N) / (PI * N).sqrt() * X.sqrt()
        for r in rs:
            out[r] += pref * (AP[r][0](X) * k0 + AP[r][1](X) * k1)
        if N == 1:
            T1 = abs(pref * (APabs[rs[0]][0](X) * k0 + APabs[rs[0]][1](X) * k1))
        tb = max(tail_bound(r, u, u, N) for r in rs)
        if tb < (T1 * arb("1e-70")).lower():
            break
        N += 1
    for r in rs:
        out[r] += arb(0, tail_bound(r, u, u, N))
    return out


def crude_sup(r, a, bnd):
    """upper bound of |D^r g| on [a, bnd]."""
    ea = fq(a).exp()
    eb = fq(bnd).exp()
    tot = arb(0)
    M = Mdeg[r]
    N = 1
    while True:
        Xa = PI * N * ea
        Xb = PI * N * eb
        k0, k1 = K01(2 * Xa)
        term = 4 * divisor_count(N) / (PI * N).sqrt() * Xb.sqrt() * (APabs[r][0](Xb) * k0 + APabs[r][1](Xb) * k1)
        tot += term
        if N == 1:
            T1 = term
        tb = tail_bound(r, a, bnd, N)
        if tb < (T1 * arb("1e-70")).lower():
            break
        N += 1
    return (tot + arb(tail_bound(r, a, bnd, N))).upper()


# ---------------- (2) g > 0 on [0, U] ----------------
def lower_bound_interval(a, bnd):
    c = (a + bnd) / 2
    rad = (bnd - a) / 2
    ra = fq(rad)
    v = point_vals(c, list(range(PT + 1)))
    g0, g1, g2 = v[0], v[1], v[2]
    # exact lower bound of the quadratic q(x) = g0 + g1 x + g2 x^2/2 on |x| <= rad
    lo_q = min((g0 + g1 * ra + g2 * ra * ra / 2).lower(), (g0 - g1 * ra + g2 * ra * ra / 2).lower())
    if g2.lower() > 0:
        x = -g1 / g2
        if not (abs(x).lower() > ra.upper()):          # vertex may lie inside
            lo_q = min(lo_q, (g0 - g1 * g1 / (2 * g2)).lower())
    elif not (g2.upper() < 0):                          # sign of g2 unknown
        lo_q = min(lo_q, (g0 - abs(g1) * ra - abs(g2) * ra * ra / 2).lower())
    hi = arb(0)
    for k in range(3, PT + 1):
        hi += abs(v[k]) * ra ** k / FACT[k]
    rem = arb(crude_sup(PT + 1, a, bnd)) * ra ** (PT + 1) / FACT[PT + 1]
    lo = arb(lo_q) - hi - rem
    return lo, g0


stack = [(Fraction(k, 64), Fraction(k + 1, 64)) for k in range(int(U * 64) - 1, -1, -1)]
nint = 0
worst = None
depth_max = 0
while stack:
    a, bnd = stack.pop()
    lo, gc = lower_bound_interval(a, bnd)
    nint += 1
    if lo.lower() > 0:
        rel = float((lo / gc).lower())
        if worst is None or rel < worst[0]:
            worst = (rel, float(a), float(bnd))
        depth_max = max(depth_max, (Fraction(1, 64) / (bnd - a)).numerator.bit_length() - 1)
        if nint % 200 == 0:
            print("    ... %d intervals, at u=%.4f, %.0fs" % (nint, float(bnd), time.time() - t0), flush=True)
        continue
    if bnd - a < Fraction(1, 2 ** 60):
        raise SystemExit("FAILED to certify g > 0 near u = %s (g(c) = %s, lower bound %s)" % (float(a), gc, lo))
    m = (a + bnd) / 2
    stack.append((m, bnd))
    stack.append((a, m))
log["intervals_[0,U]"] = nint
log["max_bisection_depth"] = depth_max
log["min_lowerbound_over_g(centre)"] = worst
print("(2) g > 0 on [0, %s] proved with %d intervals (max depth %d below 1/64, %.0fs); min (lower bound / g(centre)) = %.3g on [%.6f, %.6f]"
      % (float(U), nint, depth_max, time.time() - t0, worst[0], worst[1], worst[2]), flush=True)


# ---------------- (3) g > 0 on [U, oo) ----------------
def a_coef(nu, k):
    num = Fraction(1)
    for i in range(1, k + 1):
        num *= (4 * nu * nu - (2 * i - 1) ** 2)
    f = 1
    for i in range(1, k + 1):
        f *= i
    return num / (f * 8 ** k)


ELL = 6
P0, P1 = P[0]
M = len(P0) - 1                     # deg P0 = 4 + 2d ; deg P1 = 3 + 2d
assert M == 4 + 2 * d and len(P1) - 1 == 3 + 2 * d
XU = PI * fq(U).exp()
y = 1 / XU
theta = arb("0.5 +/- 0.5")
coef = {}
for (Pc, nu) in ((P0, 0), (P1, 1)):
    for k, pk in enumerate(Pc):
        if pk == 0:
            continue
        for i in range(ELL + 1):
            c = fq(pk) * fq(a_coef(nu, i)) / arb(2) ** i
            if i == ELL:
                c = c * theta
            e = (M - k) + i
            coef[e] = coef.get(e, arb(0)) + c
c0 = coef.pop(0)
rest = arb(0)
for e, c in coef.items():
    rest += abs(c) * y ** e
Flow = c0 - rest
ratio = arb("2.2") * Pabs1[0] * arb(2) ** M * (-2 * XU).exp() / (Flow * (1 - arb(2) ** M * (-2 * XU).exp()))
log.update(tail_leading=c0.str(10), tail_Flow=Flow.str(10), tail_ratio=ratio.str(5))
assert c0.lower() > 0 and Flow.lower() > 0 and ratio.upper() < 1, ("tail argument failed; raise U", c0, Flow, ratio)
print("(3) u >= %s (X >= %.1f): F/X^M >= %s > 0 and sum_{N>=2}|T_N| / T_1 <= %s < 1  [DLMF 10.40.10, %d terms]" %
      (float(U), float(XU.mid()), Flow.str(5), ratio.str(3), ELL), flush=True)

# ---------------- (4) cost and g(0) ----------------
g0 = point_vals(Fraction(0), [0])[0]
cost = arb(0)
NMAXC = max(80, Mdeg[0] + 2)
terms = {}
for m in range(2, NMAXC + 1):
    lam = von_mangoldt(m)
    if lam == 0:
        continue
    val = arb(0)
    N = 1
    while True:
        X = PI * N * m                       # u = log m:  X = pi N m exactly
        k0, k1 = K01(2 * X)
        pref = 4 * divisor_count(N) / (PI * N).sqrt() * X.sqrt()
        val += pref * (AP[0][0](X) * k0 + AP[0][1](X) * k1)
        if N == 1:
            T1 = abs(pref * (APabs[0][0](X) * k0 + APabs[0][1](X) * k1))
        tb = arb(0)
        N1 = max(N, Mdeg[0]) + 1
        for NN in range(N + 1, N1 + 1):
            h = arb("4.4") * Pabs1[0] * (PI * NN * m) ** Mdeg[0] * (-2 * PI * NN * m).exp()
            tb += h / (1 - RATIO) if NN == N1 else h
        tb = tb.upper()
        if tb < (T1 * arb("1e-70")).lower():
            break
        N += 1
    val += arb(0, tb)
    p_ = next(p for p in range(2, m + 1) if m % p == 0)      # m = p^k, Lambda(m) = log p
    t = 2 * arb(p_).log() / arb(m).sqrt() * val
    terms[m] = t
    cost += t
# n > NMAXC: |g(log n)| <= 4.4 Pabs(1) (pi n)^M e^{-2 pi n}/(1-0.0051)^2 ; 2 Lambda(n) n^{-1/2} <= 2 log n ;
# consecutive-n ratio of these bounds <= (1+1/n)^{M+1} e^{-2pi} <= 0.0051 for n >= M+1
Xn = PI * (NMAXC + 1)
tailc = 2 * arb(NMAXC + 1).log() * arb("4.4") * Pabs1[0] * Xn ** Mdeg[0] * (-2 * Xn).exp() / (1 - RATIO) ** 3
cost += arb(0, tailc.upper())
ratio_c = cost / g0
log.update(g0=g0.str(25), cost=cost.str(25), cost_over_g0=ratio_c.str(25),
           cost_terms={m: terms[m].str(12) for m in terms if m <= 13})
print("(4) g(0) = %s,  cost = %s" % (g0.str(20), cost.str(20)))
cup = ratio_c.upper()
print("    cost/g(0) in %s ;  certified: cost/g(0) <= %s" % (ratio_c.str(20), cup.str(12, radius=False)), flush=True)
log["certified_upper_bound"] = cup.str(15, radius=False)
log["time_s"] = time.time() - t0
json.dump(log, open("certificate_A_d%d.json" % d, "w"), indent=1)
