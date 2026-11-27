"""Attacks 3a-3c (reviewer code).  argv[1] in {a, b, c}."""
import sys; sys.path.insert(0, "/home/bas/praime/writeup/code/study/bootstrap/maass/review")
import time, json, numpy as np, weilmat as W, sdp
which = sys.argv[1]
out = []
def run(geo, r, eps, ks=(1, 2, 3, 4), tag=""):
    t0 = time.time()
    res = sdp.tstar(geo, r, eps, ks)
    print("%s N=%d X=%g eps=%d r=%.6f ks=%s t*=% .6e a_p: %s lam/k: %s %.0fs" % (
        tag, geo.N, geo.X, eps, r, "".join(map(str, ks)), res["t"],
        " ".join("% .4f" % res["ap"][p] for p in sorted(res["ap"])),
        " ".join("% .4e" % res["lam"][k] for k in ks), time.time() - t0), flush=True)
    out.append(dict(tag=tag, N=geo.N, X=geo.X, eps=eps, r=r, ks=list(ks), t=res["t"],
                    ap={str(p): v for p, v in res["ap"].items()}))
if which == "a":
    Ns = [int(v) for v in sys.argv[2].split(",")]
    for N in Ns:
        geo = W.Geo(8.0, N)
        for r in (13.0, 15.0, 18.0):
            run(geo, r, 1, tag="3a")
elif which == "b":
    geo = W.Geo(8.0, int(sys.argv[2]))
    rs = [float(v) for v in sys.argv[3].split(",")]
    for r in rs:
        run(geo, r, 0, tag="3b")
elif which == "c":
    geo = W.Geo(8.0, int(sys.argv[2]))
    for r in (9.0, 9.533695261353557, 10.0):
        for ks in ((1,), (2,), (1, 2), (1, 2, 3, 4)):
            run(geo, r, 1, ks, tag="3c")
json.dump(out, open("/home/bas/praime/writeup/code/study/bootstrap/maass/review/attack_%s_%s.json" % (which, "_".join(sys.argv[2:3])), "w"), indent=1)
