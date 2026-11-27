"""Attack 5: the tenth form (odd 19.4847) in the chain. t* near 19.48 with the key's forms and with the chain's own
fed-back forms (chain_X16_found.json, used as point masses at arccos(a^/2) as the chain does), with and without the
even neighbour 19.4235; active forms; the left crossing of -0.0005 that the rule would have used. usage: attack5.py N"""
import sys, json
sys.path.insert(0, "/home/bas/praime/writeup/code/study/bootstrap/maass/review3")
import numpy as np
from eng3 import Eng
from sdp3 import Problem
from common import load_key, bdict

N = int(sys.argv[1])
E = Eng(16.0, N)
key = load_key()
KK = [dict(r=g["r"], eps=g["eps"], b=bdict(g["ap"], E.pp)) for g in key]
ch = json.load(open("/home/bas/praime/writeup/data/study/maass/chain_X16_found.json"))
CC = [dict(r=g["r"], eps=g["eps"], b=bdict({int(p): float(np.clip(v, -2, 2)) for p, v in g["ap"].items()}, E.pp)) for g in ch]
true10 = [g for g in key if abs(g["r"] - 19.4847) < 1e-3][0]


def run(r, src, drop=None):
    kn = [k for k in src if k["r"] <= r - 0.05 + 1e-9 and (drop is None or abs(k["r"] - drop) > 0.01)]
    P = Problem(E, r, 1, True, kn)
    lb, ub = P.solve()
    x = P.x
    per = {}
    for b in P.blocks:
        if b[2] == "form":
            M = b[0] + np.tensordot(x[:-1], b[1][:-1], axes=(0, 0))
            nm = P.forms[b[3]]["name"]
            per[nm] = min(per.get(nm, np.inf), np.linalg.eigvalsh(M)[0] / b[4])
    act = [k for k, v in per.items() if v < lb + 2e-4]
    return lb, P.ap(), act, len(kn)


print("N = %d; true a_2, a_3 of 19.4847: %+.4f %+.4f" % (N, true10["ap"][2], true10["ap"][3]))
for r in [19.475, 19.480, 19.4847, 19.490]:
    row = []
    for lab, src, drop in [("key", KK, None), ("chain", CC, None), ("chain-no19.42", CC, 19.4235), ("key-no19.42", KK, 19.4235)]:
        t, ap, act, nk = run(r, src, drop)
        row.append("%s %+.5f(%d) a2 %+.4f a3 %+.4f act %s" % (lab, t, nk, ap[2], ap[3], ",".join(act)))
    print("r %.4f | " % r + " | ".join(row), flush=True)

# the rule's left crossing of t* = -0.0005 on [19.47, 19.49] with the chain's data, and the ternary search
lo, hi = 19.47, 19.49
for _ in range(10):
    mid = 0.5 * (lo + hi)
    if run(mid, CC)[0] >= -0.0005:
        hi = mid
    else:
        lo = mid
t, ap, act, nk = run(hi, CC)
print("chain data: left crossing of -0.0005 at r = %.5f (true 19.48471): a2 %+.4f a3 %+.4f a5 %+.4f a7 %+.4f (true %+.4f %+.4f %+.4f %+.4f)" % (
    hi, ap[2], ap[3], ap[5], ap[7], true10["ap"][2], true10["ap"][3], true10["ap"][5], true10["ap"][7]))
