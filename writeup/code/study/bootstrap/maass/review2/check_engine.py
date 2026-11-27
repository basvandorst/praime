"""Check weil2.Engine against zeros: Q(f) from the matrix vs sum over zeros of |f^(gamma)|^2.
f(v) = (1-(v/a)^2)^6 (in the Legendre span). Zeros from PARI (lfunzeros) and mpmath (zetazero)."""
import sys
sys.path.insert(0, "/home/bas/praime/writeup/code/study/bootstrap/maass/review2")
import numpy as np
import mpmath as mp
import cypari2
from weil2 import Engine, prime_powers, gauss, phis

pari = cypari2.Pari()
pari.allocatemem(200 * 10**6)
pari.set_real_precision(30)


def fvec(eng, deg=6):
    v, w = gauss(80, -eng.a, eng.a)
    f = (1 - (v / eng.a) ** 2) ** deg
    c = phis(eng.N, eng.a, v).T @ (w * f)
    v, w = gauss(3000, -eng.a, eng.a)   # many nodes: f^ is needed up to |t| = 160
    f = (1 - (v / eng.a) ** 2) ** deg
    return c, v, w, f


def fhat(t, v, w, f):
    return np.sum(w * f * np.exp(1j * v * t))


def check(name, X, kappas, logq, coef, zeros, pole=False):
    eng = Engine(X, 16)
    c, v, w, f = fvec(eng)
    M = eng.arch_sym(kappas) + logq * np.eye(eng.N)
    for p, m in prime_powers(X):
        cc = coef(p, m)
        if cc != 0:
            M += eng.prime_sym(p ** m, cc)
    Q = c @ M @ c
    S = sum(abs(fhat(g, v, w, f)) ** 2 for g in zeros)
    if pole:  # zeta: sum_gamma h = Q + h(i/2) + h(-i/2)
        # h(t) = f^(t) f^(-t) for real f; h(i/2) + h(-i/2) = 2 (int f e^{-v/2})(int f e^{v/2})
        hp = 2 * np.sum(w * f * np.exp(-v / 2)) * np.sum(w * f * np.exp(v / 2))
        Q = Q + hp
    print("%-28s X=%5.1f  Q(matrix) = %.14f   sum over zeros = %.14f   diff %.2e  (%d zeros)"
          % (name, X, Q, S, Q - S, len(zeros)))


def lzeros(L, T):
    z = [float(x) for x in pari.lfunzeros(L, T)]
    return [x for x in z] + [-x for x in z if x > 0]


if __name__ == "__main__":
    # zeta
    zz = [float(mp.im(mp.zetazero(n))) for n in range(1, 120)]
    zz = zz + [-x for x in zz]
    for X in (8, 16):
        check("zeta (pole)", X, [0], 0.0, lambda p, m: np.log(p), zz, pole=True)
    # chi_-3 (odd): Gamma_R(s+1), q = 3
    L3 = pari.lfuncreate(-3)
    z3 = lzeros(L3, 160)
    chi3 = lambda p: 0 if p == 3 else (1 if p % 3 == 1 else -1)
    for X in (8, 16):
        check("L(chi_-3)", X, [1], np.log(3), lambda p, m: chi3(p) ** m * np.log(p), z3)
    # chi_5 (even): Gamma_R(s), q = 5
    L5 = pari.lfuncreate(5)
    z5 = lzeros(L5, 160)
    chi5 = lambda p: 0 if p == 5 else (1 if p % 5 in (1, 4) else -1)
    for X in (8, 16):
        check("L(chi_5)", X, [0], np.log(5), lambda p, m: chi5(p) ** m * np.log(p), z5)
    # Delta, weight 12: Gamma_C(s+11/2) = Gamma_R(s+11/2) Gamma_R(s+13/2)
    mf = pari.mfinit([1, 12], 1)
    F = pari.mfbasis(mf)[0]
    LD = pari.lfunmf(mf, F)
    zD = lzeros(LD, 120)
    tau = {p: int(pari.mfcoef(F, p)) for p in range(2, 40)}

    def cD(p, m):
        lam = tau[p] / p ** 5.5
        th = np.arccos(lam / 2)
        return 2 * np.cos(m * th) * np.log(p)
    for X in (8, 16):
        check("Delta", X, [5.5, 6.5], 0.0, cD, zD)
