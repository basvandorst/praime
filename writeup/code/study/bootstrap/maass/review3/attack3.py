"""Attack 3: what sets the right ends of the pieces of round 171. For each true form j, at r = r_j + D:
 (a) the registered rule: known = true forms with r_g <= r - 0.05
 (b) known = true forms with r_g < r_j - 1e-6 (never the form j itself)
 (c) known = true forms with r_g <= r_j + 1e-6 (the form j itself always included, at distance D)
usage: attack3.py parity N"""
import sys, time
sys.path.insert(0, "/home/bas/praime/writeup/code/study/bootstrap/maass/review3")
import numpy as np
from eng3 import Eng
from sdp3 import Problem
from common import load_key, bdict

eps = int(sys.argv[1]); N = int(sys.argv[2])
key = load_key()
E = Eng(16.0, N)
KN = [dict(r=g["r"], eps=g["eps"], b=bdict(g["ap"], E.pp)) for g in key]
targets = [g for g in key if g["eps"] == eps and g["r"] > (15.5 if eps == 1 else 13.0)]
Ds = [-0.02, -0.01, 0.0, 0.01, 0.02, 0.03, 0.045, 0.06, 0.10, 0.15, 0.25]
for g in targets:
    rj = g["r"]
    for D in Ds:
        r = round(rj + D, 6)
        out = []
        for rule in "abc":
            if rule == "a":
                kn = [k for k in KN if k["r"] <= r - 0.05]
            elif rule == "b":
                kn = [k for k in KN if k["r"] < rj - 1e-6]
            else:
                if D <= 0:
                    out.append("      --      "); continue
                kn = [k for k in KN if k["r"] <= rj + 1e-6]
            P = Problem(E, r, eps, eps == 1, kn)
            lb, ub = P.solve()
            out.append("%+.5f(%d)" % (lb, len(kn)))
        print("form %.4f  D %+.3f  r %.4f   (a) %s  (b) %s  (c) %s" % (rj, D, r, *out), flush=True)
