"""Truth with the central zero: lambda_min of Q_f - |f^(0)|^2 (and Sym^2) for the odd forms, by basis size and quadrature."""
import sys, json
sys.path.insert(0, "/home/bas/praime/writeup/code/study/bootstrap/maass/review2")
import numpy as np
from sdp2 import true_value
from maass2 import maass_forms
from weil2 import Engine
key = json.load(open("/home/bas/praime/writeup/data/study/maass/lmfdb_key.json"))
X = float(sys.argv[1])
for N, nu in ((8, 800), (16, 800), (24, 800), (24, 1600), (32, 800), (40, 1200)):
    eng = Engine(X, N, nu=nu)
    row = []
    for f in key:
        if f["parity"] != "odd" or not f["ap"]:
            continue
        loc = {int(p): np.cos(np.arange(30) * np.arccos(np.clip(a / 2, -1, 1))) for p, a in f["ap"].items()}
        tc = true_value(eng, maass_forms(f["r"], 1, 2, True), loc)
        t0 = true_value(eng, maass_forms(f["r"], 1, 1, False), loc)
        row.append("%.4f: %+.2e (no centre %+.2e)" % (f["r"], tc[0], t0[0]))
    print("X=%g N=%2d nu=%d  " % (X, N, nu) + "; ".join(row), flush=True)
