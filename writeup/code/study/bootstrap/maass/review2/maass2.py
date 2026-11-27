"""Priority 3: Maass forms of SL(2,Z), parameter r, parity eps, one measure on the Satake angle per prime.
usage: maass2.py X N kmax eps center r1 r2 ...   (prints t*, the dual bound, a_2, a_3, a_5, a_7)"""
import sys, json
sys.path.insert(0, "/home/bas/praime/writeup/code/study/bootstrap/maass/review2")
import numpy as np
from sdp2 import Form, sym_weight, run, true_value
from weil2 import Engine


def maass_forms(r, eps, kmax, center):
    g = {1: [eps + 1j * r, eps - 1j * r],
         2: [0, 2j * r, -2j * r],
         3: [eps + 3j * r, eps - 3j * r, eps + 1j * r, eps - 1j * r],
         4: [0, 2j * r, -2j * r, 4j * r, -4j * r]}
    return [Form(g[k], k, sym_weight(k), center=(center and eps == 1 and k == 1)) for k in range(1, kmax + 1)]


if __name__ == "__main__":
    X, N, kmax, eps, center = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), bool(int(sys.argv[5]))
    rs = [float(v) for v in sys.argv[6:]]
    key = json.load(open("/home/bas/praime/writeup/data/study/maass/lmfdb_key.json"))
    eng = Engine(X, N)
    out = []
    for r in rs:
        forms = maass_forms(r, eps, kmax, center)
        R = run(X, N, forms, eng=eng)
        ap = {p: 2 * c[1] for p, c in R["mom"].items()}
        s = "eps=%d r=%8.4f X=%g N=%d kmax=%d center=%d  t*=%+.7f (UB %+.7f)  a_p: %s" % (
            eps, r, X, N, kmax, center, R["t"], R["UB"], " ".join("%d:%+.3f" % (p, ap[p]) for p in sorted(ap) if p <= 7))
        for f in key:
            if abs(f["r"] - r) < 5e-4 and (f["parity"] == "odd") == (eps == 1):
                loc = {int(p): np.cos(np.arange(4 * 6) * np.arccos(np.clip(a / 2, -1, 1))) for p, a in f["ap"].items()}
                tv = true_value(eng, forms, loc)
                s += "\n      truth r=%.6f: lambda_min/k per form %s; true a_p %s" % (
                    f["r"], " ".join("%+.2e" % v for v in tv), " ".join("%s:%+.3f" % (p, f["ap"][p]) for p in ("2", "3", "5", "7")))
        print(s, flush=True)
        out.append(dict(r=r, eps=eps, X=X, N=N, kmax=kmax, center=center, t=R["t"], UB=R["UB"], ap=ap))
    tag = "_".join("%g" % r for r in rs[:3])
    json.dump(out, open("/home/bas/praime/writeup/code/study/bootstrap/maass/review2/maass2_X%g_N%d_k%d_e%d_c%d_%s.json"
                        % (X, N, kmax, eps, center, tag), "w"), indent=1, default=float)
