"""Task 4: width of {r : t*(r) >= 0} around r_1 = 9.5336952613535575 at X = 8, as a function of the basis size N
(float64; valid while the peak value t*(r_1) is far above the float floor ~1e-12)."""
import sys; sys.path.insert(0, "/home/bas/praime/writeup/code/study/bootstrap/maass/review")
import time, json, numpy as np, weilmat as W, sdp
from scipy.optimize import brentq
R1 = 9.5336952613535575543
N = int(sys.argv[1]); D0 = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-11; X = float(sys.argv[2]) if len(sys.argv) > 2 else 8.0
geo = W.Geo(X, N)
log = []
def f(r):
    res = sdp.tstar(geo, r, 1)
    log.append((r, res["t"]))
    return res["t"]
t0 = time.time()
tc = f(R1)
print("N=%d X=%g t*(r1)=%.6e" % (N, X, tc), flush=True)
# expand each side until t* < 0 (t*(r1) > 0 is guaranteed by the true form)
def expand(sign):
    d = D0
    while d < 0.6:
        v = f(R1 + sign * d)
        print("  side %+d d=%.1e t*=% .6e" % (sign, d, v), flush=True)
        if v < 0:
            return d
        d *= 6.0
    raise SystemExit("no sign change within 0.6 on side %+d" % sign)
lo = expand(-1); hi = expand(+1)
# roots
rm = brentq(f, R1 - lo, R1 - lo / 6.0 if lo > 1e-11 else R1, xtol=1e-15, rtol=1e-15, maxiter=60)
rp = brentq(f, R1 + hi / 6.0 if hi > 1e-11 else R1, R1 + hi, xtol=1e-15, rtol=1e-15, maxiter=60)
print("N=%d X=%g  r_- = %.15f  r_+ = %.15f  width %.3e  centre-r1 %.3e  (%d solves, %.0fs)" % (
    N, X, rm, rp, rp - rm, 0.5 * (rm + rp) - R1, len(log), time.time() - t0), flush=True)
json.dump(dict(N=N, X=X, tc=tc, rm=rm, rp=rp, log=log), open(
    "/home/bas/praime/writeup/code/study/bootstrap/maass/review/width_X%g_N%d.json" % (X, N), "w"), indent=1)
