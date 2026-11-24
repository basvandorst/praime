"""Tabulate g(eps) = zeta'/zeta(1 + i eps) + 1/(i eps) (regular, g(0) = Euler's gamma) on the bklib grid
eps = 0 .. 7.6 step 0.0025 with mpmath (25 digits). usage: gtab.py -> gtab.npz"""
import os, numpy as np, mpmath as mp
mp.mp.dps = 25
eps = np.arange(0, 7.6 + 0.00125, 0.0025); g = np.empty(len(eps), complex)
for i, e in enumerate(eps):
    if e == 0: g[i] = float(mp.euler); continue
    s = mp.mpc(1, e); v = mp.zeta(s, derivative=1) / mp.zeta(s) + 1 / (1j * mp.mpf(e)); g[i] = complex(v)
np.savez(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'gtab.npz'), eps=eps, g=g)
# check against a series: zeta'/zeta(1+x) = -1/x + gamma + (-gamma^2 - 2 gamma_1) x + ...
x = 1j * 1e-3; ser = complex(mp.euler) + (-float(mp.euler) ** 2 - 2 * float(mp.stieltjes(1))) * x
print('g(0.001) table %.12f%+.12fj  series %.12f%+.12fj' % (g[0].real if False else complex(mp.zeta(1 + x, derivative=1) / mp.zeta(1 + x) + 1 / x).real, complex(mp.zeta(1 + x, derivative=1) / mp.zeta(1 + x) + 1 / x).imag, ser.real, ser.imag))
print('n =', len(eps), 'g(1) =', g[400], 'g(7.6) =', g[-1])
