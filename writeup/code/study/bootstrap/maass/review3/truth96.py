"""Truth at large basis: each true form f (key r, key a_p as point masses) against L(f) (central zero if odd), Sym^2 and
L(f x g) for every true g with r_g <= r_f - 0.05 (rule of round 171). Smallest eigenvalue / k of each form, worst one."""
import sys
sys.path.insert(0, "/home/bas/praime/writeup/code/study/bootstrap/maass/review3")
import numpy as np
from eng3 import Eng
from sdp3 import true_min
from common import load_key, bdict

key = load_key()
for N in [int(a) for a in sys.argv[1:]]:
    E = Eng(16.0, N)
    KK = [dict(r=g["r"], eps=g["eps"], b=bdict(g["ap"], E.pp)) for g in key]
    for f in key:
        kn = [k for k in KK if k["r"] <= f["r"] - 0.05 + 1e-9]
        tm = true_min(E, f["r"], f["eps"], f["eps"] == 1, f["ap"], kn)
        worst = min(tm, key=tm.get)
        print("N %d form %.4f %s: %d RS forms; worst %s %+.2e; L %+.2e S2 %+.2e; RS min %+.2e" % (
            N, f["r"], "odd " if f["eps"] else "even", len(kn), worst, tm[worst], tm["L"], tm["S2"],
            min([v for k, v in tm.items() if k.startswith("RS")] or [np.nan])), flush=True)
