# kernels k(z) = z^(-1/2) K_3(z): Meijer G (mpmath) against the Mellin-Barnes integral of the gamma ratio (numpy), three signatures;
# and the degree-2 conventions against the Bessel kernels of the cubic case.
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, mpmath as mp, math, pred4 as P
from scipy.special import j0, y0, k0
mp.mp.dps = 20
def k_meijer(z, mus):
    a = [(mp.mpf(1) / 2 + m) / 2 for m in mus]; d = len(mus)
    return float(2 * mp.meijerg([[], []], [a, [1 - x for x in a]], (mp.pi**d * z)**2)) / math.sqrt(z)
def k_mb(z, mus, c=0.1, T=20000.0, dt=0.004):
    t = np.arange(0, T, dt); w = c + 1j * t; g = P.Kt(w, mus) * np.exp(-w * math.log(z))
    wts = np.ones(len(t)); wts[0] = wts[-1] = 0.5
    return float((g * wts).sum().real * dt / math.pi)
worst = 0.0
for r2, mus in P.SIG_MU.items():
    for z in (0.02, 0.3, 1.0, 3.7, 12.0, 40.0):
        a, b = k_meijer(z, mus), k_mb(z, mus); worst = max(worst, abs(a - b))
        print('r2=%d mu=%s z=%6.2f  meijer %+.9f  MB %+.9f  diff %.1e' % (r2, mus, z, a, b, abs(a - b)))
for z in (0.05, 0.7, 2.3):
    u = 4 * math.pi * math.sqrt(z)
    print('deg 2 D<0: meijer %+.9f  2pi J0 %+.9f | D>0: meijer %+.9f  2pi(-Y0+2K0/pi) %+.9f' % (k_meijer(z, (0, 1)), 2 * math.pi * j0(u), k_meijer(z, (0, 0)), 2 * math.pi * (-y0(u) + 2 / math.pi * k0(u))))
print('largest |meijer - MB| (degree 3): %.1e' % worst)
