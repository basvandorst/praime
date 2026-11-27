"""Attack 6: the even piece at 19.4235. Round 171 has {19.43} only; eng3/sdp3 at N = 32 find t* > 0 up to 19.47.
Basis convergence and active forms. usage: attack6.py eps N r1,r2,..."""
import sys
sys.path.insert(0, "/home/bas/praime/writeup/code/study/bootstrap/maass/review3")
import numpy as np
from eng3 import Eng
from sdp3 import Problem
from common import load_key, bdict

eps = int(sys.argv[1]); N = int(sys.argv[2]); rs = [float(v) for v in sys.argv[3].split(",")]
excl = float(sys.argv[4]) if len(sys.argv) > 4 else None   # rule (b): never use the form at excl
E = Eng(16.0, N)
key = load_key()
KK = [dict(r=g["r"], eps=g["eps"], b=bdict(g["ap"], E.pp)) for g in key]
for r in rs:
    kn = [k for k in KK if k["r"] <= r - 0.05 + 1e-9 and (excl is None or abs(k["r"] - excl) > 1e-3)]
    P = Problem(E, r, eps, eps == 1, kn)
    lb, ub = P.solve()
    x = P.x
    per = {}
    for b in P.blocks:
        if b[2] == "form":
            M = b[0] + np.tensordot(x[:-1], b[1][:-1], axes=(0, 0))
            nm = P.forms[b[3]]["name"]
            per[nm] = min(per.get(nm, np.inf), np.linalg.eigvalsh(M)[0] / b[4])
    ap = P.ap()
    print("eps %d N %d r %.4f excl %s known %d: t* in [%+.6f, %+.6f]; a_p %s; per form %s" % (
        eps, N, r, excl, len(kn), lb, ub, " ".join("%d:%+.3f" % (p, ap[p]) for p in sorted(ap)),
        " ".join("%s:%+.5f" % kv for kv in per.items())), flush=True)
