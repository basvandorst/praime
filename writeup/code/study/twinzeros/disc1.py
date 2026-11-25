"""disc1.py -- discovery range: pipeline checks on theorems (explicit formula), then the first pattern fits."""
import numpy as np, sys, json
import tzlib, theory

ZDIR = "/home/bas/praime/writeup/data/riemann/"
k, X, cells = tzlib.load_cells(["disc.npz"])
cum = tzlib.cumulative(cells)
Xf = X.astype(float)
u = np.log(Xf)
sq = np.sqrt(Xf)
j0 = int(np.searchsorted(X, 10 ** 6))
j1 = int(np.searchsorted(X, 2 * 10 ** 9, side="right")) - 1
print("window u = %.4f .. %.4f, %d grid points" % (u[j0], u[j1], j1 - j0 + 1))


def comb(gammas, uu, chunk=200):
    """-2 Re sum_{gamma>0} e^{i gamma u}/rho for a real character (zeros symmetric), rho = 1/2 + i gamma."""
    out = np.zeros(len(uu))
    for i in range(0, len(gammas), chunk):
        g = gammas[i:i + chunk]
        rho = 0.5 + 1j * g
        out += -2 * np.real(np.exp(1j * np.outer(uu, g)) @ (1 / rho))
    return out


def load_zeros(name, tmax):
    z = np.loadtxt(ZDIR + name) if name.endswith(".txt") else None
    z = z[(z > 0) & (z <= tmax)]
    return z


Tc = 1500.0
zz = load_zeros("zeros_40000_flint.txt", Tc)
z3 = load_zeros("lzeros/D-3.txt", Tc)
z4 = load_zeros("lzeros/D-4.txt", Tc)
z5 = load_zeros("lzeros/D5.txt", Tc)
print("zeros up to %.0f: zeta %d, chi_-3 %d, chi_-4 %d, chi_5 %d" % (Tc, len(zz), len(z3), len(z4), len(z5)))
uu = u.copy()
Czeta = comb(zz, uu); C3 = comb(z3, uu); C4 = comb(z4, uu); C5 = comb(z5, uu)

chars = theory.primitive_characters(30)
byl = {c["label"]: c for c in chars}
def Ychi(label, doubled=False):
    c = byl[label]
    return tzlib.char_series(cum, c["v"], c["q"], doubled) / sq

Yth = (cum["th"] - Xf) / sq
Y3 = Ychi("3.1").real
Y4 = Ychi("4.1").real
Y5 = Ychi("5.1").real
# class 2 mod 3: theta(x;3,2) - x/2
Y32 = (cum["c3"][:, 2] - Xf / 2) / sq

GLO, GHI = 2.0, 1000.0
print("\n-- theorem checks (explicit formula), band [%g, %g], zeros to %g --" % (GLO, GHI, Tc))
for name, y, regs, pred in (("theta on zeta comb", Yth, [Czeta], [1]),
                             ("theta(chi_-3) on its comb", Y3, [C3], [1]),
                             ("theta(chi_-4) on its comb", Y4, [C4], [1]),
                             ("theta(chi_5) on its comb", Y5, [C5], [1]),
                             ("theta(x;3,2)-x/2 on zeta, chi_-3 combs", Y32, [Czeta, C3], [0.5, -0.5])):
    b, cov, d = tzlib.fit(y, regs, u, j0, j1, GLO, GHI, 1)
    print("  %-40s beta = %s  err = %s  pred %s  chi2/dof %.2f" % (name, np.round(b, 4), np.round(np.sqrt(np.diag(cov)), 4), pred, d["chi2_dof"]))

print("\n-- twins h=2 on the regressors from primes, band [%g, %g] --" % (GLO, GHI))
S2 = theory.singular_const("h2")
YP = (cum["P2"] - S2 * Xf) / sq
for regs_names in (["zeta"], ["zeta", "3.1"], ["zeta", "3.1", "4.1", "5.1"]):
    regs = []
    for r in regs_names:
        regs.append(Yth if r == "zeta" else Ychi(r).real)
    b, cov, d = tzlib.fit(YP, regs, u, j0, j1, GLO, GHI, 2)
    print("  regs %-28s beta = %s err = %s chi2/dof %.2f  (pred zeta 2S = %.4f, 3.1: 0, 4.1: 0, 5.1: 2S/3 = %.4f)" % (
        regs_names, np.round(b, 3), np.round(np.sqrt(np.diag(cov)), 3), d["chi2_dof"], 2 * S2, 2 * S2 / 3))
for regs_names, combs in ((["zeta comb"], [Czeta]), (["zeta", "chi-3", "chi-4", "chi5 combs"], [Czeta, C3, C4, C5])):
    b, cov, d = tzlib.fit(YP, combs, u, j0, j1, GLO, GHI, 2)
    print("  combs %-28s beta = %s err = %s chi2/dof %.2f" % (regs_names, np.round(b, 3), np.round(np.sqrt(np.diag(cov)), 3), d["chi2_dof"]))
