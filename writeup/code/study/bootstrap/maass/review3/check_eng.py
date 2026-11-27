"""Checks of eng3: (1) explicit formula of zeta with its zeros; (2) A_kappa against the t-domain integral
(1/pi) int h(t) Re psi_R(1/2 + kappa + it) dt for complex kappa (Maass and Rankin-Selberg sizes)."""
import sys
sys.path.insert(0, "/home/bas/praime/writeup/code/study/bootstrap/maass/review3")
import numpy as np
import mpmath as mp
from scipy.special import digamma
from eng3 import Eng, gl, basis, EULER

X, N = 16.0, 24
E = Eng(X, N, nu=700)
a = E.a
vq, wq = gl(400, -a, a)


def coeffs(fun):
    B = basis(N, a, vq)
    return (B * (wq * fun(vq))[:, None]).sum(0)


def fhat(fun, t):
    t = np.atleast_1d(t)
    return (np.exp(1j * np.outer(t, vq)) * (wq * fun(vq))[None]).sum(1)


rng = np.random.default_rng(1)
P = rng.normal(size=4)
funs = {"bump4": lambda v: (1 - (v / a) ** 2) ** 4,
        "bump4*poly": lambda v: (1 - (v / a) ** 2) ** 4 * (P[0] + P[1] * v / a + P[2] * (v / a) ** 2 + P[3] * (v / a) ** 3),
        "bump5 odd": lambda v: (1 - (v / a) ** 2) ** 5 * (v / a)}

# (1) zeta
mp.mp.dps = 20
NZ = 120
gam = np.array([float(mp.zetazero(n).imag) for n in range(1, NZ + 1)])
Mz = E.arch([0.0])
for p, m in E.pp:
    Mz = Mz + E.prime(p, m)          # c(p^m) = log p (the log p is inside prime())
for name, fun in funs.items():
    c = coeffs(fun)
    lhs = c @ Mz @ c
    hz = np.abs(fhat(fun, gam)) ** 2
    Fp = (wq * fun(vq) * np.exp(-vq / 2)).sum(); Fm = (wq * fun(vq) * np.exp(vq / 2)).sum()
    rhs = 2 * hz.sum() - 2 * Fp * Fm
    print("zeta  %-12s  matrix %.12f   zeros-poles %.12f   diff %.2e  (last zero term %.1e)" % (name, lhs, rhs, lhs - rhs, 2 * hz[-1]))

# (2) archimedean term in the t-domain
tt, wt = [], []
for lo in np.arange(-400, 400, 2.0):
    x, w = gl(40, lo, lo + 2.0)
    tt.append(x); wt.append(w)
tt = np.concatenate(tt); wt = np.concatenate(wt)


def psiR(s):
    return -0.5 * np.log(np.pi) + 0.5 * digamma(s / 2)


for kap in [0.0, 1.0, 1 + 16.1381j, 2j * 16.1381, 1 + 30.5j, 0.35j, 1 + 0.06j]:
    M = E.arch([kap, np.conj(kap)]) if np.imag(kap) != 0 else E.arch([kap])
    for name, fun in funs.items():
        c = coeffs(fun)
        lhs = c @ M @ c
        h = np.abs(fhat(fun, tt)) ** 2
        ks = [kap, np.conj(kap)] if np.imag(kap) != 0 else [kap]
        rhs = sum((wt * h * np.real(psiR(0.5 + k + 1j * tt))).sum() / np.pi for k in ks)
        print("arch kappa=%-14s %-12s matrix %.12f  t-domain %.12f  diff %.2e" % (kap, name, lhs, rhs, lhs - rhs))

# (3) convergence in nu
for nu in (350, 1400):
    E2 = Eng(X, N, nu=nu)
    d = np.max(np.abs(E2.arch([1 + 35.0j, 1 - 35.0j]) - E.arch([1 + 35.0j, 1 - 35.0j])))
    print("nu %d vs 700: max |diff| of arch(1 +- 35i) = %.2e" % (nu, d))
