"""disc3.py -- discovery range: all patterns, full model (zeta, all predicted characters jointly, null characters)."""
import numpy as np, sys
import tzlib, theory, models

D = models.Data(["disc.npz"])
GLO, GHI = 2.0, float(sys.argv[1]) if len(sys.argv) > 1 else 6000.0
per_cond = len(sys.argv) > 2 and sys.argv[2] == "percond"
for name in ("h2", "h4", "h6", "h8", "h10", "h12", "h30", "SG", "T026", "T046"):
    top = 2 * 10 ** 9 if name != "SG" else 2 * 10 ** 9
    j0, j1 = D.J(10 ** 6), D.J(top, right=True) - 1
    Y, S = D.Ypat(name)
    S_, regs = models.build(D, name, per_conductor=per_cond)
    b, cov, d = tzlib.fit(Y, [r[1] for r in regs], D.u, j0, j1, GLO, GHI, models.kforms(name))
    e = np.sqrt(np.diag(cov))
    print("\n== %s  S = %.5f  band [%g, %g]  window 1e6..%.0e  chi2/dof %.3f" % (name, S, GLO, GHI, top, d["chi2_dof"]))
    nulls = []
    for (lab, ser, pred), x, s in zip(regs, b, e):
        z = (x - pred) / s
        print("   %-28s measured %8.3f +- %6.3f   predicted %7.3f   z = %+5.2f" % (lab, x, s, pred, z))
        if lab.startswith("null"):
            nulls.append(z)
    if nulls:
        # joint chi^2 of the null coefficients with their covariance
        idx = [i for i, r in enumerate(regs) if r[0].startswith("null")]
        bb = b[idx]; C = cov[np.ix_(idx, idx)]
        chi2 = float(bb @ np.linalg.solve(C, bb))
        print("   null set: chi2 = %.2f on %d dof" % (chi2, len(idx)))
