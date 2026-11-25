"""Transform of murmuration densities against the family-theorem kernel.

(a) Bessel-Mellin identity (GR 6.561.14) at s = -2it:  int_0^inf x^{-2it} J_nu(x) dx = 2^{-2it} Gamma((nu+1-2it)/2)/Gamma((nu+1+2it)/2),
    checked with mpmath.quadosc for nu = 1, 3 and several t.
(b) E(M) := -(1/2pi) int_0^inf y^{-1} M(y) G(y) dy for M = M_harm, M_zub, M_prime and Gaussian H (centre t0, width sig),
    against the predicted line sums:
      harmonic : -(1/2pi) sum_m mu^2(m)/phi(m) Hhat(2 ln m)                       (main term: int H)
      Zubrilina: -(1/2pi) sum_{d,s} (alpha/2pi) Q(d) d/s^2 Hhat(2 ln(d/s))           (main term: gamma/(2 pi^2) int H)
      prime nat: -(1/2pi) sum_{(s,d)=1} A((s/d)^2) Hhat(2 ln(d/s)), A = (6/pi^2) s^{-2} prod_{p|d} p^2/((p-1)^2(p+1)) (main: 6/pi^2)
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import mpmath as mp
import sympy as sp
from densities import M_harm, M_zub, M_prime, G_kernel, H_int, H_hat, ALPHA, GAMMA

# (a)
mp.mp.dps = 20
print("(a) Bessel-Mellin identity int_0^inf x^{-2it} J_nu(x) dx:")
for nu in (1, 3):
    for t in (0.3, 2.0, 7.5):
        f = lambda x: x ** (-2j * t) * mp.besselj(nu, x)
        num = mp.quadosc(f, [0, mp.inf], omega=1)
        ex = mp.mpf(2) ** (-2j * t) * mp.gamma((nu + 1 - 2j * t) / 2) / mp.gamma((nu + 1 + 2j * t) / 2)
        print(f"   nu={nu} t={t}: numeric {mp.nstr(num, 12)}  exact {mp.nstr(ex, 12)}  |diff| {mp.nstr(abs(num - ex), 3)}")

# (b)
def predicted_lines(kind, t0, sig, smax=60):
    tot = 0j
    if kind == "harm":
        for m in range(1, 2000):
            f = sp.factorint(m)
            if any(e > 1 for e in f.values()):
                continue
            phi = 1
            for q in f:
                phi *= q - 1
            tot += H_hat(2 * np.log(m), t0, sig) / phi
    elif kind == "zub":
        for d in range(1, 400):
            f = sp.factorint(d)
            if any(e > 1 for e in f.values()):
                continue
            Q = 1.0
            for q in f:
                Q *= q ** 2 / (q ** 4 - 2 * q ** 2 - q + 1)
            for s in range(1, smax):
                tot += ALPHA / (2 * np.pi) * Q * d / s ** 2 * H_hat(2 * np.log(d / s), t0, sig)
    elif kind == "prime":
        for d in range(1, 400):
            f = sp.factorint(d)
            if any(e > 1 for e in f.values()):
                continue
            w = 1.0
            for q in f:
                w *= q ** 2 / ((q - 1) ** 2 * (q + 1))
            for s in range(1, smax):
                if np.gcd(s, d) != 1:
                    continue
                tot += 6 / np.pi ** 2 * w / s ** 2 * H_hat(2 * np.log(d / s), t0, sig)
    return -tot / (2 * np.pi)


for (t0, sig) in ((10.0, 2.0), (14.0, 2.0), (18.0, 2.5)):
    ymax = ((t0 + 7 * sig) / (2 * np.pi)) ** 2 * 1.6
    ymin = max(1e-4, ((t0 - 7 * sig) / (2 * np.pi)) ** 2 / 3) if t0 - 7 * sig > 0 else 1e-4
    y = np.linspace(ymin, ymax, 60001)
    G = G_kernel(y, t0, sig)
    dy = y[1] - y[0]
    norm = -H_int(t0, sig) / (2 * np.pi)
    print(f"\n(b) Gaussian H: t0={t0}, sigma={sig}, y in [{ymin:.4f}, {ymax:.2f}]; values / (-int H/2pi):")
    for name, M, kind, main in (("harmonic", M_harm, "harm", 1.0), ("Zubrilina k=2", M_zub, "zub", GAMMA / (2 * np.pi ** 2)),
                                ("prime natural (derived)", M_prime, "prime", 6 / np.pi ** 2)):
        Mv = M(y)
        E = -(1 / (2 * np.pi)) * np.sum(Mv * G / y) * dy
        pl = predicted_lines(kind, t0, sig)
        print(f"   {name:26s}: transform {E / norm:.6f}   predicted line sum {pl / norm:.6f}   main term only {main:.6f}")
