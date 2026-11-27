"""Attack 2/3: t* for the odd problem with the central zero, Sym^2 and L(f x g) for the four lowest true forms
(9.5337, 12.1730, 13.7798, 14.3585; key data), exact moments, own barrier solver, several basis sizes."""
import sys, time
sys.path.insert(0, "/home/bas/praime/writeup/code/study/bootstrap/maass/review3")
import numpy as np
from eng3 import Eng
from sdp3 import Problem, true_min
from common import load_key, bdict

key = load_key()
low4 = [k for k in key if k["r"] < 15]
rs = [16.00, 16.10, 16.1381, 16.13807317152103, 16.20, 16.30]
Ns = [int(a) for a in sys.argv[1:]] or [24, 32, 48]
for N in Ns:
    E = Eng(16.0, N)
    known = [dict(r=g["r"], eps=g["eps"], b=bdict(g["ap"], E.pp)) for g in low4]
    for r in rs:
        t0 = time.time()
        P = Problem(E, r, 1, True, known)
        lb, ub = P.solve()
        ap = P.ap()
        print("N %2d r %.5f  t* in [%+.6f, %+.6f]  a2 %+.4f a3 %+.4f a5 %+.4f a7 %+.4f  (toep min %.1e, %.1f s)" % (
            N, r, lb, ub, ap[2], ap[3], ap[5], ap[7], P.toep, time.time() - t0), flush=True)
    # the true fifth form against the same problem (point masses at its angles)
    f5 = [k for k in key if abs(k["r"] - 16.1381) < 1e-3][0]
    tm = true_min(E, f5["r"], 1, True, f5["ap"], known)
    print("N %2d  true form 16.1381, smallest eigenvalue / k of each form:" % N,
          " ".join("%s %+.2e" % (k, v) for k, v in tm.items()), flush=True)
