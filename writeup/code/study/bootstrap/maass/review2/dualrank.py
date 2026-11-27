"""Rank of the dual matrices (test-function mixtures) at the optimum, per parity block; and a sign-free decoupling:
each prime power p^m its own measure with weight x (not x^m), coefficient in [-1,1] log p."""
import sys
sys.path.insert(0, "/home/bas/praime/writeup/code/study/bootstrap/maass/review2")
import numpy as np
from sdp2 import Form, char_weight, sym_weight, build, solve, certify, run
from maass2 import maass_forms
from weil2 import Engine


def ranks(eng, forms, label):
    var, blocks, meta = build(eng, forms)
    x, Z, gap = solve(var, blocks, meta)
    sc = sum(meta[i][1].k * np.trace(Z[i]) for i in range(len(blocks)) if meta[i][0] == "form")
    out = []
    for i in range(len(blocks)):
        if meta[i][0] == "form":
            ev = np.sort(np.linalg.eigvalsh(Z[i] / sc))[::-1]
            out.append("k=%d %s: %s" % (meta[i][1].k, "even" if meta[i][2][0] == 0 else "odd",
                                         " ".join("%.1e" % v for v in ev[:4])))
    print(label, "t*=%.9f gap %.1e |" % (x[-1], gap), " | ".join(out), flush=True)


ranks(Engine(14, 32), [Form([1], 1, char_weight)], "chi odd X=14 N=32")
ranks(Engine(24, 32), [Form([0], 1, char_weight)], "chi even X=24 N=32")
e16 = Engine(16, 40)
ranks(e16, maass_forms(9.533695261353557, 1, 2, True), "Maass odd 9.5337 centre X=16 N=40")
ranks(e16, maass_forms(14.358509518, 1, 2, True), "Maass odd 14.3585 centre X=16 N=40")

# sign-free decoupling: each (p, m) one measure, coefficient E[x] log p in [-1, 1] log p
signfree = lambda p, m: {1: 1.0}
for X in (8, 14, 24):
    for kap, q in ((1, 3), (0, 5)):
        R = run(X, 32, [Form([kap], 1, signfree)], decouple=True)
        R2 = run(X, 32, [Form([kap], 1, char_weight)], decouple=True)
        print("X=%d kappa=%d: smallest conductor, one measure per prime power with weight x^m: %.6f (cert %.6f); "
              "with weight x (sign free at even m): %.6f (cert %.6f)" % (X, kap, np.exp(-R2["t"]), np.exp(-R2["UB"]),
                                                                       np.exp(-R["t"]), np.exp(-R["UB"])), flush=True)
