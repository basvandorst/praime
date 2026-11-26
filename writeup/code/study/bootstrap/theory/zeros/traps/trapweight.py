"""D3, weight part (rigorous): recompute (ii) and the weight lower bound of trapcert_k<k>.json with a sharp upper
bound of h on [a_k, b_k] (trapcert.py used h(c) + |h'(c)| r + sup|h''| r^2/2 with the crude coefficientwise sup of
|h''|, which is useless for the wide outer traps because h has coefficients up to 1e22).
Here: [a_k, b_k] split into 2^m pieces; on each, Taylor polynomial of order 16 at the centre with ball derivatives,
exact maximum of the quadratic part, |terms 3..16| <= |D^j h(c)| r^j/j!, remainder by crude_sup(D^17 h) r^17/17!
(the mirror image of certify.py's lower bound).  Same h, eta', a_k, b_k as the certificate.
weight >= (eta W(h) - kt (1 + (eta h(0) - 1)^+)) / (eta sup h),  eta = eta'/g0,  kt = W(gt)/g0  (balls).
Usage: python trapweight.py k [m=4]   -> trapweight_k<k>.json"""
import sys, json, os, time
from fractions import Fraction
import common
from flint import arb, ctx
ctx.prec = 800
from fcert import Cert, fq   # noqa: E402

k = int(sys.argv[1])
m = int(sys.argv[2]) if len(sys.argv) > 2 else 4
t0 = time.time()
T = json.load(open("trapcert_k%d.json" % k))
C = json.load(open(os.path.join(common.CERTS, "certificate_A_d24.json")))
bg = [Fraction(s) for s in C["b_rational"]]
bh = [Fraction(s) for s in T["h_rational"]]
etaP = Fraction(T["eta_prime"])
a_k, b_k = Fraction(T["a_k"]), Fraction(T["b_k"])
CH = Cert(bh, PT=16)


def upper_bound_interval(Cx, a, bnd):
    PT = Cx.PT
    c = (a + bnd) / 2
    ra = fq((bnd - a) / 2)
    v = Cx.point_vals(c, list(range(PT + 1)))
    g0, g1, g2 = v[0], v[1], v[2]
    hi_q = max((g0 + g1 * ra + g2 * ra * ra / 2).upper(), (g0 - g1 * ra + g2 * ra * ra / 2).upper())
    if g2.upper() < 0:
        x = -g1 / g2
        if not (abs(x).lower() > ra.upper()):
            hi_q = max(hi_q, (g0 - g1 * g1 / (2 * g2)).upper())
    elif not (g2.lower() > 0):
        hi_q = max(hi_q, (g0 + abs(g1) * ra + abs(g2) * ra * ra / 2).upper())
    rest = arb(0)
    for j in range(3, PT + 1):
        rest += abs(v[j]) * ra ** j / Cx.FACT[j]
    rem = arb(Cx.crude_sup(PT + 1, a, bnd)) * ra ** (PT + 1) / Cx.FACT[PT + 1]
    return (arb(hi_q) + rest + rem).upper()


pieces = [(a_k + (b_k - a_k) * i / 2 ** m, a_k + (b_k - a_k) * (i + 1) / 2 ** m) for i in range(2 ** m)]
ups = [upper_bound_interval(CH, lo, hi) for (lo, hi) in pieces]
suph = max(ups)
CG = Cert(bg, PT=1)
g0 = CG.value0()
kt = CG.cost_ball() / g0
h0 = CH.value0()
Wh = CH.cost_ball()
eta = fq(etaP) / g0
x = eta * h0 - 1
pos = arb(0) if x.upper() < 0 else (x if x.lower() > 0 else arb(0).union(x.upper()))
num = eta * Wh - kt * (1 + pos)
assert num.lower() > 0, "(ii) fails"
supw = eta * arb(suph)
wlb = num.lower() / supw.upper()
out = dict(k=k, pieces=2 ** m, sup_h_upper=arb(suph).str(15, radius=False), ii_margin=num.str(15),
           weight_lower=arb(wlb).str(12, radius=False), W_h=Wh.str(30), h0=h0.str(20), kt=kt.str(25), g0=g0.str(25),
           eta=eta.str(20), time_s=time.time() - t0)
json.dump(out, open("trapweight_k%d.json" % k, "w"), indent=1)
print("k=%d: sup_[a_k,b_k] h <= %s ; eta W(h) - kt(...) = %s ; weight >= %s (%.0fs)" % (
    k, arb(suph).str(12, radius=False), num.str(10), arb(wlb).str(10, radius=False), time.time() - t0))
