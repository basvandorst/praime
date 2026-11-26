# Checks of kerlib.py: (1) Mellin transforms of the kernels vs the gamma ratios; (2) Mellin-Barnes inversion for K_+;
# (3) Poisson-summed primitives C(Z) vs the direct Bessel series (smooth cutoff); (4) large-Z asymptotics vs exact.
import math, numpy as np, mpmath as mp
from scipy.special import j1, y1, k1
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kerlib as K
mp.mp.dps = 12
Jk = lambda u: mp.besselj(0, u); Kp = lambda u: -mp.bessely(0, u) + 2 / mp.pi * mp.besselk(0, u)
MJ = lambda s: 2**s * mp.gamma((1 + s) / 2) / mp.gamma((1 - s) / 2)
MK = lambda s: 4**s * mp.gamma((1 + s) / 4)**2 / mp.gamma((1 - s) / 4)**2
print('(1)-(2) Mellin-Barnes: int_0^oo K(u) u^a exp(-u^2) du (scipy)  vs  (1/2pi i) int_(c) M(s) Gamma((a-s)/2)/2 ds (mpmath)')
from scipy.integrate import quad
from scipy.special import j0, y0, k0
KJ = lambda u: j0(u); KP = lambda u: -y0(u) + 2 / math.pi * k0(u)
for a in (0.5, 1.0, 2.0, 3.0):
    for name, Kf, M in (('J0', KJ, MJ), ('K+', KP, MK)):
        direct = quad(lambda u: Kf(u) * u**a * math.exp(-u * u), 0, 12, limit=400, points=[1e-6, 1e-3, 0.1, 1])[0]
        vals = []
        for c in (-0.5, 0.2):
            f = lambda t: mp.re(M(c + 1j * t) * mp.gamma((a - c - 1j * t) / 2) / 2)
            vals.append(float(mp.quad(f, [-mp.inf, -20, 0, 20, mp.inf]) / (2 * mp.pi)))
        print('  a=%.1f %s: direct %.12f   MB(c=-0.5) %.12f   MB(c=0.2) %.12f' % (a, name, direct, vals[0], vals[1]))
print('(3) C(Z): Poisson closed form vs direct sum_n sqrt(Z) k1(4 pi n sqrt Z)/(2 pi n), smooth cutoff exp(-(n/N)^2)')
for Z in (0.05, 0.2, 0.6, 1.7, 5.3, 13.1):
    x = 4 * math.pi * math.sqrt(Z)
    for N in (20000, 80000):
        n = np.arange(1, 6 * N, dtype=float); w = np.exp(-(n / N)**2)
        dJ = np.sum(math.sqrt(Z) * j1(n * x) / (2 * math.pi * n) * w)
        dK = np.sum(math.sqrt(Z) * (-y1(n * x) - 2 / math.pi * k1(n * x)) / (2 * math.pi * n) * w)
    print('  Z=%5.2f  C-: direct %.6f  closed %.6f   C+: direct %.6f  closed %.6f' % (Z, dJ, K.C_minus_exact([Z])[0], dK, K.C_plus_exact([Z])[0]))
print('    small-Z forms at Z=0.05, 0.2: C-', K.C_small([0.05, 0.2], -1), ' C+', K.C_small([0.05, 0.2], 1))
print('(4) large Z: exact vs asymptotic')
Zs = np.array([10.3, 25.7, 59.9, 60.1, 150.2, 400.9])
for sgn in (-1, 1):
    ex = (K.C_minus_exact if sgn < 0 else K.C_plus_exact)(Zs); asy = K.C_asym(Zs, sgn)
    print('  sgn %+d  max |exact - asym| = %.2e  (values %s)' % (sgn, np.max(np.abs(ex - asy)), np.round(ex, 4)))
    for Z in (59.0, 61.0, 99.3):
        print('     Z=%.1f exact %.6f asym %.6f' % (Z, (K.C_minus_exact if sgn < 0 else K.C_plus_exact)(np.array([Z]))[0], K.C_asym(np.array([Z]), sgn)[0]))
