"""Task 2: t*(r) at X = 8, kmax = 4, exact moment relaxation, Legendre basis N (argv[1])."""
import sys; sys.path.insert(0, "/home/bas/praime/writeup/code/study/bootstrap/maass/review")
import time, json, numpy as np, weilmat as W, sdp
N = int(sys.argv[1]); X = float(sys.argv[2]) if len(sys.argv) > 2 else 8.0
pts = [(1, r) for r in (3.0, 6.0, 9.0, 9.53369526, 10.0, 11.0, 12.5)] + [(0, r) for r in (12.0, 13.77975135, 15.0)]
geo = W.Geo(X, N)
out = []
for eps, r in pts:
    t0 = time.time()
    res = sdp.tstar(geo, r, eps)
    rec = dict(X=X, N=N, eps=eps, r=r, t=res["t"], gap=res["gap"], ap={str(p): v for p, v in res["ap"].items()},
               lam={str(k): v for k, v in res["lam"].items()}, toepmin=res["toepmin"])
    out.append(rec)
    print("eps=%d r=%-11s t*=% .6f (gap %.0e) a_p: %s  per-form lam/k: %s  %.0fs" % (
        eps, r, res["t"], res["gap"], " ".join("% .4f" % res["ap"][p] for p in sorted(res["ap"])),
        " ".join("% .5f" % res["lam"][k] for k in (1, 2, 3, 4)), time.time() - t0), flush=True)
json.dump(out, open("/home/bas/praime/writeup/code/study/bootstrap/maass/review/task2_X%g_N%d.json" % (X, N), "w"), indent=1)
