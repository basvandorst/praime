"""Truth check with the forced central zero removed: true odd forms (a_p for p <= 31 from the key; X = 8 needs p <= 7),
lambda_min of (M_f - h(0) e0 e0^T), and of M_f; X = 8 and 16."""
import sys; sys.path.insert(0, "/home/bas/praime/writeup/code/study/bootstrap/maass/review")
import json, numpy as np, weilmat as W
key = json.load(open("/home/bas/praime/writeup/data/study/maass/lmfdb_key.json"))
for X in (8.0, 16.0):
    for N in (24, 48):
        geo = W.Geo(X, N)
        for f in key:
            if f["parity"] != "odd" or f["r"] > 19: continue
            ap = {int(p): v for p, v in f["ap"].items()}
            if max(p for p, m in W.prime_powers(X)) > max(ap): continue
            mom = W.true_moments(ap)
            M = W.form_matrix(geo, 1, f["r"], 1, mom)
            Mc = M.copy(); Mc[0, 0] -= 2 * geo.a
            print("X=%g N=%d r=%.4f  lam_min M_f = % .3e   with central zero removed = % .3e" % (
                X, N, f["r"], np.linalg.eigvalsh(M)[0], np.linalg.eigvalsh(Mc)[0]), flush=True)
