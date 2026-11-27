"""Normalisation checks of weilmat.py (reviewer code).
(i) zeta: full explicit formula with the first 150 zeros (mpmath.zetazero) vs the x-space matrix.
(ii) t-space vs x-space for complex shifts mu (L(f) odd r=9.5337, Sym^4 r=9.5337), smooth phi.
"""
import numpy as np, mpmath as mp
import sys; sys.path.insert(0, "/home/bas/praime/writeup/code/study/bootstrap/maass/review")
import weilmat as W

mp.mp.dps = 20
X, N = 8.0, 16
geo = W.Geo(X, N)
a = geo.a
# smooth test function phi(x) = (1-(x/a)^2)^3, coefficients in the orthonormal Legendre basis
xq, wq = W.gl(60, -a, a)
phi = (1 - (xq / a) ** 2) ** 3
B = W.basis_vals(N, a, xq)
c = B.T @ (wq * phi)
print("projection residual", np.max(np.abs(B @ c - phi)))


def phihat(t):
    # a int_{-1}^1 (1-s^2)^3 cos(a t s) ds = a sqrt(pi) Gamma(4) (2/z)^(7/2) J_{7/2}(z), z = a t
    z = mp.mpf(a) * t
    if abs(z) < 1e-3:
        return float(np.sum(wq * phi * np.cos(float(t) * xq)))
    return float(mp.mpf(a) * mp.sqrt(mp.pi) * 6 * (2 / z) ** 3.5 * mp.besselj(3.5, z))


# (i) zeta
mz = W.arch(geo, [0])
for p, m in W.prime_powers(X):
    mz -= 2 * np.log(p) * p ** (-m / 2) * geo.G(m * np.log(p))
Qx = c @ mz @ c
pole = 2 * float(np.sum(wq * phi * np.cosh(xq / 2))) ** 2
zs = 0.0
for n in range(1, 151):
    g = float(mp.im(mp.zetazero(n)))
    zs += 2 * phihat(g) ** 2
print("pole", pole, "form", Qx, "check phihat(3)", phihat(3.0), float(np.sum(wq * phi * np.cos(3.0 * xq))))
print("zeta: zero side %.12f   x-space form + pole %.12f   diff %.2e" % (zs, Qx + pole, zs - Qx - pole))

# (ii) t-space with complex shifts


def tspace(mus):
    f = lambda t: phihat(t) ** 2 * sum(2 * mp.re(W.psiR(mp.mpf(0.5) + mp.mpc(mu) + 1j * t)) for mu in mus)
    # even integrand: (1/2pi) int_R = (1/pi) int_0^inf
    pts = [0] + [float(k) for k in range(5, 200, 5)] + [mp.inf]
    return float(mp.quad(f, pts) / mp.pi)


for name, mus in [("L(f) eps=1 r=9.5337", W.gammas(1, 9.53369526, 1)),
                  ("L(f) eps=0 r=13.78", W.gammas(1, 13.77975135, 0)),
                  ("Sym2 r=9.5337", W.gammas(2, 9.53369526, 1)),
                  ("Sym4 r=9.5337", W.gammas(4, 9.53369526, 1))]:
    Ax = c @ W.arch(geo, mus) @ c
    At = tspace(mus)
    print("%-22s x-space %.12f  t-space %.12f  diff %.2e" % (name, Ax, At, Ax - At))
