"""Priority 2: level 1, holomorphic weight kappa as a real number, one measure on the Satake angle per prime.
Forms L(f), Sym^2, Sym^3, Sym^4 (Gamma factors as derived in review2.md), conductor e^{-kt} for Sym^k.
usage: hol2.py X N kmax center(0/1) kappa1 kappa2 ..."""
import sys, json
sys.path.insert(0, "/home/bas/praime/writeup/code/study/bootstrap/maass/review2")
import numpy as np
from sdp2 import Form, sym_weight, run, true_value
from weil2 import Engine, primes_upto


def hol_forms(kap, kmax, center):
    g = {1: [(kap - 1) / 2, (kap + 1) / 2],
         2: [1.0, kap - 1, kap],
         3: [3 * (kap - 1) / 2, 3 * (kap - 1) / 2 + 1, (kap - 1) / 2, (kap + 1) / 2],
         4: [0.0, kap - 1, kap, 2 * kap - 2, 2 * kap - 1]}
    return [Form(g[k], k, sym_weight(k), center=(center and k == 1)) for k in range(1, kmax + 1)]


import cypari2
_pari = cypari2.Pari()
TRUE = {}
for _k in (12, 16, 18, 20, 22):
    _F = _pari.mfbasis(_pari.mfinit([1, _k], 1))[0]
    TRUE[_k] = {p: int(_pari.mfcoef(_F, p)) for p in range(2, 40) if all(p % q for q in range(2, p))}

X, N, kmax, center = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), bool(int(sys.argv[4]))
kaps = [float(v) for v in sys.argv[5:]]
eng = Engine(X, N)
out = []
for kap in kaps:
    forms = hol_forms(kap, kmax, center)
    R = run(X, N, forms, eng=eng)
    ap = {p: 2 * c[1] for p, c in R["mom"].items()}
    line = dict(kappa=kap, X=X, N=N, kmax=kmax, center=center, t=R["t"], UB=R["UB"], ap=ap)
    s = "kappa=%6.2f X=%g N=%d kmax=%d center=%d  t*=%+.7f (UB %+.7f)  a_p: %s" % (
        kap, X, N, kmax, center, R["t"], R["UB"], " ".join("%d:%+.4f" % (p, a) for p, a in ap.items()))
    ki = int(round(kap))
    if abs(kap - ki) < 1e-12 and ki in TRUE:
        lam = {p: TRUE[ki][p] / p ** ((ki - 1) / 2) for p in primes_upto(X)}
        loc = {p: np.cos(np.arange(4 * 6) * np.arccos(lam[p] / 2)) for p in lam}
        tv = true_value(eng, forms, loc)
        line["true_per_form"] = tv
        line["true_ap"] = lam
        s += "\n        truth: per-form lambda_min/k = %s; true a_p: %s" % (
            " ".join("%+.2e" % v for v in tv), " ".join("%d:%+.4f" % (p, a) for p, a in lam.items()))
    print(s, flush=True)
    out.append(line)
json.dump(out, open("/home/bas/praime/writeup/code/study/bootstrap/maass/review2/hol2_X%g_N%d_k%d_c%d.json"
                    % (X, N, kmax, center), "w"), indent=1, default=float)
