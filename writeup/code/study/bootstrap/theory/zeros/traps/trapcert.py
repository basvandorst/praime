"""D3: rigorous trap certificate at degree 48 (half-degree d' = 24) for node k (ball arithmetic, 800 bits).

gt  = the certified cushioned function of certs/certificate_A_d24.json (b_rational, exact rationals).  ASSUMED from
      that certificate (certify.py, cert_d24.out): B_gt > 0 on R and gt > 0 on [0, oo), i.e. gt is in the cone.
      gt(0) = g0 is a ball near 1; the lemma is applied to gt/g0 (cost kt = W(gt)/g0, a ball computed here).
h   = h_k of hp_d24.json rounded to SIGH significant digits, top coefficient set to 0, minus DELTA * g_0
      (g_0 = Phi*Phi > 0: lowers h everywhere, helps every condition except (ii), which keeps a wide margin).
eta'= rational (FAC * eta_0 rounded to 3 digits);  the lemma's eta is eta'/g0.  With F := gt - eta' h:
  (i)   eta B_h <= B_{gt/g0}  <=>  B_F = B_gt - eta' B_h >= 0 : exact Sturm count of q_F(s) on [0, oo) + q_F(0) > 0;
  (iii) eta b_24(h) = 0 <= b_24(gt)/g0 : exact;
  S     w = eta h - gt/g0 = -F/g0, so S = {F < 0} ; we prove F > 0 on [0, a_k] and [b_k, U] (Taylor models,
        bisection, as certify.py) and on [U, oo) (asymptotic argument of certify.py), hence S c (a_k, b_k);
  (ii)  eta W(h) - kt (1 + (eta h(0) - 1)^+) > 0 with balls (W by the cost code of certify.py);
  weight bound: sup_S w <= eta sup_{[a_k, b_k]} h (gt >= 0), sup h <= h(c) + |h'(c)| r + sup|h''| r^2/2.
Usage: python trapcert.py k [FAC=2] [DELTA=1/1000] [SIGH=90] [U=6]   -> trapcert_k<k>.json, trapcert_k<k>.log"""
import sys, json, time, os, math
from decimal import Decimal, getcontext
from fractions import Fraction
import common
from flint import arb, fmpq, ctx
import sympy

ctx.prec = 800
from fcert import Cert, fq          # noqa: E402
import fam                          # noqa: E402  (numerics for the endpoints only)

k = int(sys.argv[1])
FAC = Fraction(sys.argv[2]) if len(sys.argv) > 2 else Fraction(2)
DELTA = Fraction(sys.argv[3]) if len(sys.argv) > 3 else Fraction(1, 1000)
SIGH = int(sys.argv[4]) if len(sys.argv) > 4 else 90
U = Fraction(sys.argv[5]) if len(sys.argv) > 5 else Fraction(6)
tag = "k%d" % k + ("" if FAC == 2 else "_x%s" % FAC)
LOG = open("trapcert_%s.log" % tag, "w")


def log(s):
    print(s, flush=True)
    LOG.write(s + "\n"); LOG.flush()


t0 = time.time()
ctx.prec = 800
C = json.load(open(os.path.join(common.CERTS, "certificate_A_d24.json")))
bg = [Fraction(s) for s in C["b_rational"]]
assert len(bg) == 25 and bg[24] > 0
R = json.load(open("hp_d24.json"))
H = R["h"][k - 1]
getcontext().prec = SIGH


def rnd(s):
    x = +Decimal(s)                    # rounded to SIGH significant digits
    return Fraction(x)


bh = [rnd(s) for s in H["b"]]
bh[24] = Fraction(0)
bh[0] -= DELTA
uk_s = R["u"][k - 1]
lam_f, W_f = float(R["lam"]), float(H["W"])
eta0_f = lam_f / W_f
e = FAC * Fraction(eta0_f)
ex = math.floor(math.log10(float(e)))
etaP = Fraction(round(float(e) / 10 ** ex, 2)).limit_denominator(1000) * Fraction(10) ** ex
nk = min([q for q in range(2, 40) if len([p for p in range(2, q + 1) if q % p == 0 and all(p % r for r in range(2, p))]) == 1],
         key=lambda q: abs(math.log(q) - float(uk_s)))
log("node k=%d (n_k=%d), eta' = %s (= %.4f eta_0), DELTA = %s, SIGH = %d, U = %s" % (k, nk, etaP, float(etaP) / eta0_f, DELTA, SIGH, U))
bF = [bg[j] - etaP * bh[j] for j in range(25)]

# ---------- numerics: endpoints of {F < 0} around u_k ----------
ctx.prec = 700
Fn = fam.Fam(24, 1)
bFa = [fq(x) for x in bF]
PF0 = Fn.comb(bFa, 0)
PF1 = Fn.comb(bFa, 1)
uk = arb(uk_s)
gpp = arb(R["gpp"][k - 1])


def Fv(u, r=0):
    return Fn.eval_polys([PF0 if r == 0 else PF1], u.exp())[0]


assert Fv(uk) < 0
r0 = (2 * fq(etaP) / gpp).sqrt()
ends = []
for s in (-1, 1):
    L = r0 * arb("1.5")
    while Fv(uk + s * L) < 0:
        L = L * 2
    lo, hi = arb(0), L
    x = r0
    for it in range(100):
        ux = fam.mid(uk + s * x)
        f0 = Fv(ux)
        f1 = Fv(ux, 1) * s
        if f0 < 0:
            lo = x
        else:
            hi = x
        xn = fam.mid(x - f0 / f1)
        if not (xn > lo and xn < hi):
            xn = fam.mid((lo + hi) / 2)
        if abs(float(xn - x)) < 1e-30 * float(r0):
            x = xn
            break
        x = xn
    ends.append(x)
dig = int(-math.log10(float(r0))) + 15


def tofrac(x):
    return Fraction(round(Fraction(x.mid().str(dig + 30, radius=False)) * 10 ** dig), 10 ** dig)


WID = Fraction(11, 10)
a_k = tofrac(uk - ends[0] * fq(WID))
b_k = tofrac(uk + ends[1] * fq(WID))
log("numerical S around u_k: u_k - %s .. u_k + %s ; certified box (a_k, b_k) = S widened by 10%%" %
    (ends[0].str(6, radius=False), ends[1].str(6, radius=False)))
ctx.prec = 800

out = dict(k=k, n_k=nk, eta_prime=str(etaP), eta_prime_over_eta0=float(etaP) / eta0_f, DELTA=str(DELTA), SIGH=SIGH, U=str(U),
           a_k=str(a_k), b_k=str(b_k), h_rational=[str(x) for x in bh])

# ---------- (i) exact Sturm ----------
s = sympy.symbols("s")
q = sympy.Poly([sympy.Rational(x.numerator, x.denominator) for x in reversed(bF)], s, domain="QQ")
nroots = q.count_roots(0, None)
assert nroots == 0 and bF[0] > 0, ("(i) fails", nroots)
log("(i) q_F(s) = q_gt(s) - eta' q_h(s) has %d roots in [0,oo) and q_F(0) > 0: eta B_h < B_(gt/g0) on R [exact]" % nroots)
out["i_sturm_roots"] = int(nroots)
# ---------- (iii) ----------
assert bh[24] == 0 and bg[24] > 0
log("(iii) b_24(h) = 0 <= b_24(gt) [exact]")

# ---------- F > 0 on [0, a_k] u [b_k, U] and on [U, oo) ----------
CF = Cert(bF, PT=16)
log("    polynomials of D^r F built (%.0fs)" % (time.time() - t0))
ints = []
for i in range(int(U * 64)):
    lo, hi = Fraction(i, 64), Fraction(i + 1, 64)
    if hi <= a_k or lo >= b_k:
        ints.append((lo, hi))
    else:
        if lo < a_k:
            ints.append((lo, a_k))
        if hi > b_k:
            ints.append((b_k, hi))
res = CF.positive_on(ints, log=log, t0=t0)
log("(S) F > 0 on [0, a_k] and [b_k, %s]: %d intervals, max depth %d, worst lower/F(centre) = %.3g on [%.9f, %.9f] (%.0fs)" %
    (U, res["intervals"], res["depth_max"], res["worst"][0], res["worst"][1], res["worst"][2], time.time() - t0))
ok, tl = CF.tail_positive(U)
assert ok, ("tail failed", tl)
log("(S) F > 0 on [%s, oo): leading %s, Flow %s, N>=2 ratio %s" % (U, tl["tail_leading"], tl["tail_Flow"], tl["tail_ratio"]))
out.update(S_intervals=res["intervals"], S_depth=res["depth_max"], S_worst=res["worst"], tail=tl)

# ---------- (ii) and the weight bound ----------
CG = Cert(bg, PT=1)
g0 = CG.value0()
Wg = CG.cost_ball()
kt = Wg / g0
CH = Cert(bh, PT=1)
h0 = CH.value0()
Wh = CH.cost_ball()
eta = fq(etaP) / g0
x = eta * h0 - 1
if x.upper() < 0:
    pos = arb(0)
elif x.lower() > 0:
    pos = x
else:
    pos = arb(0).union(x.upper())
num = eta * Wh - kt * (1 + pos)
assert num.lower() > 0, ("(ii) fails", num)
log("(ii) g0 = %s, kt = W(gt)/g0 = %s, W(h) = %s, h(0) = %s, eta = %s" % (g0.str(15), kt.str(15), Wh.str(15), h0.str(10), eta.str(10)))
log("     eta W(h) - kt (1 + (eta h(0) - 1)^+) = %s > 0" % num.str(10))
c = (a_k + b_k) / 2
r = (b_k - a_k) / 2
hv = CH.point_vals(c, [0, 1])
s2 = CH.crude_sup(2, a_k, b_k)
suph = hv[0] + abs(hv[1]) * fq(r) + arb(s2) * fq(r) ** 2 / 2
supw = eta * suph
wlb = num.lower() / supw.upper()
logn = arb(nk).log()
Rk = max(abs(fq(a_k) - logn).upper(), abs(fq(b_k) - logn).upper())
log("weight: sup_S w <= eta sup_[a_k,b_k] h = %s ;  mu(S) >= %s  (mu_k = %s)" % (supw.str(10), arb(wlb).str(10, radius=False), R["mu"][k - 1][:14]))
log("RESULT k=%d n=%d: S c (a_k, b_k), a_k - log n = %s, b_k - log n = %s, radius R_k <= %s, weight >= %s" %
    (k, nk, (fq(a_k) - logn).str(5), (fq(b_k) - logn).str(5), arb(Rk).str(4, radius=False), arb(wlb).str(8, radius=False)))
out.update(g0=g0.str(30), kt=kt.str(30), W_h=Wh.str(30), h0=h0.str(20), eta=eta.str(20), ii_margin=num.str(15),
           sup_h=suph.str(15), sup_w_upper=supw.upper().str(10, radius=False),
           a_minus_logn=(fq(a_k) - logn).str(10), b_minus_logn=(fq(b_k) - logn).str(10),
           R_k_upper=arb(Rk).str(6, radius=False), weight_lower=arb(wlb).str(12, radius=False),
           mu_k=R["mu"][k - 1][:20], time_s=time.time() - t0,
           assumed=["gt in the cone (certify.py on certificate_A_d24.json: Sturm for B_gt, gt > 0 on [0,oo))",
                    "closed Bessel form of Phi*Phi (README of certs)", "DLMF 10.40.10 remainder bound",
                    "kbessel.K01 series ratio bounds", "python-flint/arb and sympy count_roots correct",
                    "the trap lemma itself (specification)"])
json.dump(out, open("trapcert_%s.json" % tag, "w"), indent=1)
log("done %.0fs" % (time.time() - t0))
