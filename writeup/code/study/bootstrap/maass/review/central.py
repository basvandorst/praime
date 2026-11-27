"""Not registered: odd forms have root number -1, so L(1/2, f) = 0.  Odd problem with the forced central zero removed
(Q_f - h(0) >= 0).  Compare t* with and without, X = 8."""
import sys; sys.path.insert(0, "/home/bas/praime/writeup/code/study/bootstrap/maass/review")
import time, json, numpy as np, weilmat as W, sdp
N = int(sys.argv[1]); X = float(sys.argv[2]); rs = [float(v) for v in sys.argv[3].split(",")]
geo = W.Geo(X, N); out = []
for r in rs:
    t0 = time.time(); res = sdp.tstar(geo, r, 1, central=True)
    out.append((r, res["t"], res["ap"]))
    print("central-zero odd X=%g N=%d r=%.6f t*=% .6e a_p: %s %.0fs" % (X, N, r, res["t"], " ".join("% .4f" % res["ap"][p] for p in sorted(res["ap"])), time.time() - t0), flush=True)
json.dump([(r, t, {str(p): v for p, v in ap.items()}) for r, t, ap in out], open("/home/bas/praime/writeup/code/study/bootstrap/maass/review/central_X%g_N%d_%s.json" % (X, N, sys.argv[4] if len(sys.argv) > 4 else "a"), "w"))
