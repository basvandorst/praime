import json
import sys
sys.path.insert(0, "/home/bas/praime/writeup/code/study/bootstrap/maass/review3")
from sdp3 import cheb_b

KEY = "/home/bas/praime/writeup/data/study/maass/lmfdb_key.json"


def load_key():
    d = json.load(open(KEY))
    out = []
    for e in d:
        if not e["ap"]:
            continue
        out.append(dict(r=float(e["r_str"]), eps=1 if e["parity"] == "odd" else 0,
                        ap={int(p): float(v) for p, v in e["ap_str"].items()}))
    return out


def bdict(ap, pp, eps_g=None):
    """b_m(g) = alpha^m + alpha^-m from a_p (point mass at the Satake angle)."""
    return {(p, m): cheb_b(ap[p], m)[m] for p, m in pp}
