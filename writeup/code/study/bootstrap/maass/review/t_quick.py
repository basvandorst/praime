import sys; sys.path.insert(0, "/home/bas/praime/writeup/code/study/bootstrap/maass/review")
import time, numpy as np, weilmat as W, sdp
for N in (16, 32):
    geo = W.Geo(8.0, N)
    for r in (3.0, 9.53369526):
        t0 = time.time()
        res = sdp.tstar(geo, r, 1, verbose=False)
        print(N, r, "t*=%.6f gap %.1e" % (res["t"], res["gap"]), {p: round(v, 4) for p, v in res["ap"].items()},
              {k: "%.6f" % v for k, v in res["lam"].items()}, "toep %.1e" % res["toepmin"], "%.1fs" % (time.time() - t0), flush=True)
