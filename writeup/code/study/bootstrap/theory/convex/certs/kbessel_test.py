"""Check kbessel.K01 against mpmath (150 digits) on x in [6, 3000]: the mpmath value must lie in the ball."""
from flint import arb, ctx, fmpq
import mpmath as mp
from kbessel import K01
ctx.prec = 320
mp.mp.dps = 150
print("const_euler" in dir(arb))
bad = 0; worst = 0; cnt = 0
import random
random.seed(1)
xs = [fmpq(random.randint(600, 300000), 100) for _ in range(800)] + [fmpq(k, 10) for k in range(590, 610)]
for xq in xs:
    x = arb(xq)
    k0, k1 = K01(x)
    for nu, k in ((0, k0), (1, k1)):
        ref = mp.besselk(nu, mp.mpf(int(xq.p)) / int(xq.q))
        refa = arb(mp.nstr(ref, 145))
        cnt += 1
        if not k.overlaps(refa):
            bad += 1
            print("BAD", float(x.mid()), nu, k.str(15), mp.nstr(ref, 15))
        rel = float((k.rad() / abs(k.mid())).mid()) if k.mid() != 0 else 1
        worst = max(worst, rel)
print("tested", cnt, "bad", bad, "worst relative radius", worst)
