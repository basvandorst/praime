"""Kernels.  K_d(z) = (1/2pi) int G(t) z^{it} dt,  G(t) = prod_j Gamma_R(1/2-it+mu_j)/Gamma_R(1/2+it+mu_j)
(inverse Mellin transform on Re w = 0 of the gamma ratio of Q^{it} X(1/2+it)); as a Meijer G-function
K_d(z) = 2 G^{d,0}_{0,2d}((pi^d z)^2 | a_1..a_d ; 1-a_1..1-a_d),  a_j = (1/2 + mu_j)/2.
GL(2) weight k: mu = ((k-1)/2, (k+1)/2)   -> should equal 2 pi sqrt(z) J_{k-1}(4 pi sqrt z).
Sym^2 of weight k: mu = (1, k-1, k).
K21(z) = 4 pi i^{-k} sqrt(z) cos(4 pi z) J_{k-1}(4 pi z)   (Petersson c = N^2 term)."""
import numpy as np
import mpmath as mp
from scipy.special import jv


def Kd_meijer(z, mus):
    a = [(0.5 + m) / 2 for m in mus]
    d = len(mus)
    x = (mp.pi ** d * z) ** 2
    return 2 * mp.meijerg([[], []], [a, [1 - aj for aj in a]], x)


def K2(z, k):
    return 2 * np.pi * np.sqrt(z) * jv(k - 1, 4 * np.pi * np.sqrt(z))


def K3(z, k):
    return float(Kd_meijer(z, [1, k - 1, k]))


def K21(z, k):
    z = np.asarray(z, dtype=float)
    return 4 * np.pi * (-1) ** (k // 2) * np.sqrt(z) * np.cos(4 * np.pi * z) * jv(k - 1, 4 * np.pi * z)


if __name__ == "__main__":
    mp.mp.dps = 30
    for k in (2, 4):
        print("GL(2) check, k=%d:" % k, [("%.6f" % float(Kd_meijer(z, [(k - 1) / 2, (k + 1) / 2])), "%.6f" % K2(z, k)) for z in (0.1, 0.7, 2.3)])
    for k in (2, 4):
        print("k=%d  z: K3  K21" % k)
        for z in (0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8, 16):
            print("   %6.2f %+.5f %+.5f" % (z, K3(z, k), K21(z, k)))
