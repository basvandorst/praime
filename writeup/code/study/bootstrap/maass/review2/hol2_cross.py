"""Smallest admissible weight: bisection for t*(kappa) = 0 near 12 (no central zero) and 18 (central zero).
usage: hol2_cross.py X N kmax"""
import sys
sys.path.insert(0, "/home/bas/praime/writeup/code/study/bootstrap/maass/review2")
import numpy as np
from sdp2 import Form, sym_weight, run
from weil2 import Engine


def hol_forms(kap, kmax, center):
    g = {1: [(kap - 1) / 2, (kap + 1) / 2],
         2: [1.0, kap - 1, kap],
         3: [3 * (kap - 1) / 2, 3 * (kap - 1) / 2 + 1, (kap - 1) / 2, (kap + 1) / 2],
         4: [0.0, kap - 1, kap, 2 * kap - 2, 2 * kap - 1]}
    return [Form(g[k], k, sym_weight(k), center=(center and k == 1)) for k in range(1, kmax + 1)]


X, N, kmax = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
eng = Engine(X, N)
for center, lo, hi in ((False, 11.9, 12.05), (True, 17.9, 18.05)):
    f = lambda k: run(X, N, hol_forms(k, kmax, center), eng=eng)
    Rlo, Rhi = f(lo), f(hi)
    print("center=%d X=%g N=%d kmax=%d: t*(%.4f)=%+.3e  t*(%.4f)=%+.3e" % (center, X, N, kmax, lo, Rlo["t"], hi, Rhi["t"]), flush=True)
    if Rlo["t"] * Rhi["t"] > 0:
        continue
    for it in range(30):
        mid = 0.5 * (lo + hi)
        R = f(mid)
        if R["t"] < 0:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-7:
            break
    ap = {p: round(2 * c[1], 4) for p, c in R["mom"].items()}
    k0 = round(0.5 * (lo + hi))
    Rk = f(float(k0))
    print("   smallest admissible weight %.7f (t* at %d = %+.3e, UB %+.3e); a_p at the crossing %s; at %d: %s"
          % (0.5 * (lo + hi), k0, Rk["t"], Rk["UB"], ap, k0,
             {p: round(2 * c[1], 4) for p, c in Rk["mom"].items()}), flush=True)
