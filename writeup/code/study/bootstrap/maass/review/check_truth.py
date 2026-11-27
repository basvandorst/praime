"""Task 1: smallest eigenvalues of the four Weil forms at X = 8 for true Maass forms (true r, true a_p as point masses).
Also: wrong parity and shifted r with the same a_p (power of the check)."""
import sys; sys.path.insert(0, "/home/bas/praime/writeup/code/study/bootstrap/maass/review")
import json, numpy as np
import weilmat as W

key = json.load(open("/home/bas/praime/writeup/data/study/maass/lmfdb_key.json"))
forms = {1: key[0], 3: key[2]}
X = 8.0
for N in (8, 16, 24, 32, 48, 64):
    geo = W.Geo(X, N)
    for num, f in forms.items():
        eps = 1 if f["parity"] == "odd" else 0
        mom = W.true_moments({int(p): v for p, v in f["ap"].items()})
        row = []
        for k in (1, 2, 3, 4):
            ev = np.linalg.eigvalsh(W.form_matrix(geo, k, f["r"], eps, mom))
            row.append(ev[0])
        wrong = []
        for k in (1, 3):
            wrong.append(np.linalg.eigvalsh(W.form_matrix(geo, k, f["r"], 1 - eps, mom))[0])
        shift = [np.linalg.eigvalsh(W.form_matrix(geo, 1, f["r"] + dr, eps, mom))[0] for dr in (-0.5, 0.5)]
        print("N=%2d form %d r=%.6f eps=%d  lam_min L,S2,S3,S4: %s | wrong eps L,S3: %s | L at r-+0.5: %s"
              % (N, num, f["r"], eps, " ".join("%.3e" % v for v in row), " ".join("%.3e" % v for v in wrong),
                 " ".join("%.3e" % v for v in shift)), flush=True)
