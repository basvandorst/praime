"""Priority 1: real characters, one measure on [-1,1] per prime, weight <x^m> log p at p^m. Exact moments."""
import sys, json
sys.path.insert(0, "/home/bas/praime/writeup/code/study/bootstrap/maass/review2")
import numpy as np
from sdp2 import Form, char_weight, run, true_value, atoms
from weil2 import Engine, primes_upto

chi = {1: (3, lambda p: 0 if p == 3 else (1 if p % 3 == 1 else -1)),
       0: (5, lambda p: 0 if p == 5 else (1 if p % 5 in (1, 4) else -1))}

cases = [(8, 16), (8, 32), (14, 16), (14, 32), (24, 24), (24, 40), (40, 32)]
if len(sys.argv) > 1:
    cases = [tuple(int(v) for v in a.split(",")) for a in sys.argv[1:]]
res = []
for kap in (1, 0):
    q, ch = chi[kap]
    for X, N in cases:
        eng = Engine(X, N)
        forms = [Form([kap], 1, char_weight)]
        R = run(X, N, forms, eng=eng)
        loc = {p: np.cos(np.arange(8) * np.arccos(ch(p))) for p in primes_upto(X)}
        tv = true_value(eng, forms, loc)[0]
        Rd = run(X, N, forms, decouple=True, eng=eng)
        line = dict(kappa=kap, X=X, N=N, t=R["t"], LB=R["LB"], UB=R["UB"], gap=R["gap"],
                    q_star=float(np.exp(-R["t"])), q_star_UB=float(np.exp(-R["UB"])),
                    true_lmin=tv, q_true=float(np.exp(-tv)), delta=tv + np.log(q),
                    t_minus_true=R["t"] - tv, UB_minus_true=R["UB"] - tv,
                    q_decoupled=float(np.exp(-Rd["t"])), q_dec_UB=float(np.exp(-Rd["UB"])))
        at = {}
        for p, c in R["mom"].items():
            th, ev = atoms(c)
            at[int(p)] = dict(m1=float(c[1]), Ex=None, atoms_x=None if th is None else [float(np.cos(v)) for v in th])
        # power moments E x^m
        for p, c in R["mom"].items():
            from math import comb
            Ex = [sum(comb(m, j) * c[abs(m - 2 * j)] for j in range(m + 1)) / 2 ** m for m in range(1, len(c))]
            at[int(p)]["Ex"] = [round(float(e), 6) for e in Ex]
        line["local"] = at
        line["dual_argmax_x"] = {str(k): round(float(np.cos(v)), 6) for k, v in R["argmax"].items()}
        res.append(line)
        print("kappa=%d X=%2d N=%2d  q*=%.12f (UB-cert %.12f)  true: lmin=%.12f q_true=%.12f delta=%.2e  "
              "t*-true=%.2e UB-true=%.2e  decoupled q*=%.6f (cert %.6f)"
              % (kap, X, N, line["q_star"], line["q_star_UB"], tv, line["q_true"], line["delta"],
                 line["t_minus_true"], line["UB_minus_true"], line["q_decoupled"], line["q_dec_UB"]), flush=True)
        print("    E[x^m] per prime:", {p: v["Ex"] for p, v in at.items()})
        print("    atoms (x):", {p: (None if v["atoms_x"] is None else [round(a, 4) for a in v["atoms_x"]]) for p, v in at.items()})
        print("    dual argmax x:", line["dual_argmax_x"], flush=True)
json.dump(res, open("/home/bas/praime/writeup/code/study/bootstrap/maass/review2/chars2_%s.json"
                    % "_".join("%d-%d" % c for c in cases), "w"), indent=1)
