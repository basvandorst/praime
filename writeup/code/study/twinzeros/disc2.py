"""disc2.py -- twins h=2: dependence on the band and the window; empirical null for the error bars."""
import numpy as np
import tzlib, theory

k, X, cells = tzlib.load_cells(["disc.npz"])
cum = tzlib.cumulative(cells)
Xf = X.astype(float); u = np.log(Xf); sq = np.sqrt(Xf)
chars = theory.primitive_characters(30)
byl = {c["label"]: c for c in chars}
Ychi = lambda lab: tzlib.char_series(cum, byl[lab]["v"], byl[lab]["q"]) / sq
Yth = (cum["th"] - Xf) / sq
S2 = theory.singular_const("h2")
YP = (cum["P2"] - S2 * Xf) / sq
Y3, Y4, Y5 = Ychi("3.1").real, Ychi("4.1").real, Ychi("5.1").real
J = lambda x: int(np.searchsorted(X, x))
j0, j1 = J(10 ** 6), J(2 * 10 ** 9) - 1

print("band dependence (window 1e6..2e9); columns: zeta (pred %.3f), 3.1 (0), 4.1 (0), 5.1 (pred %.3f)" % (2 * S2, 2 * S2 / 3))
for glo, ghi in ((2, 50), (2, 100), (2, 200), (2, 500), (2, 1000), (2, 2000), (2, 4000), (2, 6000),
                 (2, 200), (200, 500), (500, 1000), (1000, 2000), (2000, 4000), (4000, 6000)):
    b, cov, d = tzlib.fit(YP, [Yth, Y3, Y4, Y5], u, j0, j1, glo, ghi, 2)
    e = np.sqrt(np.diag(cov))
    print("  [%5g,%5g] bins %5d  " % (glo, ghi, d["nbins"]) + "  ".join("%6.3f+-%5.3f" % (x, y) for x, y in zip(b, e)))

print("\nwindow dependence (band [2, 1000])")
for xa, xb in ((1e5, 2e9), (1e6, 2e9), (1e7, 2e9), (1e8, 2e9), (1e6, 1e8), (1e8, 2e9), (1e4, 1e6)):
    a, bb = J(int(xa)), J(int(xb)) - 1
    b, cov, d = tzlib.fit(YP, [Yth, Y3, Y4, Y5], u, a, bb, 2, 1000, 2)
    e = np.sqrt(np.diag(cov))
    print("  [%.0e, %.0e]  " % (xa, xb) + "  ".join("%6.3f+-%5.3f" % (x, y) for x, y in zip(b, e)))

print("\nempirical null: regress the twins on the prime regressors shifted in u by s cells (s = 1000..8000 step 350)")
bs = []
for s in range(1000, 8001, 350):
    sh = lambda a: np.roll(a, s)
    # shifted regressors must stay inside the computed range: use j0..j1 with regressor taken at j - s
    regs = [np.concatenate([np.zeros(s), r[:-s]]) for r in (Yth, Y3, Y4, Y5)]
    if j0 - s < 0:
        continue
    b, cov, d = tzlib.fit(YP, regs, u, j0, j1, 2, 1000, 2)
    bs.append((b, np.sqrt(np.diag(cov))))
B = np.array([x[0] for x in bs]); E = np.array([x[1] for x in bs])
print("  %d shifts; mean beta %s ; sd of beta %s ; mean reported err %s" % (len(B), np.round(B.mean(0), 3), np.round(B.std(0, ddof=1), 3), np.round(E.mean(0), 3)))
print("  z-scores sd (should be 1):", np.round((B / E).std(0, ddof=1), 2))
