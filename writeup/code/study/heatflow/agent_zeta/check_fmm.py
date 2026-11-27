import numpy as np, time, sys
from scipy.integrate import quad
from flow import load, Force, tail_force

M = int(sys.argv[1]) if len(sys.argv) > 1 else 12000
x, Y0 = load(M)
F = Force(M, Y0)
t0 = time.time(); f = F.pair_sum(x); t1 = time.time()
print('M=%d Y0=%.3f  fmm time %.3f s' % (M, Y0, t1 - t0))
# direct, in chunks
d = np.empty(M)
for a in range(0, M, 500):
    xa = x[a:a + 500, None]
    m = 1.0 / (xa - x[None, :])
    m[np.arange(len(xa)), np.arange(a, a + len(xa))] = 0.0
    d[a:a + 500] = m.sum(1) + (1.0 / (xa + x[None, :])).sum(1)
print('direct time %.2f s; max |fmm-direct| = %.3e ; rms %.3e ; max |direct| %.3f'
      % (time.time() - t1, np.max(np.abs(f - d)), np.sqrt(np.mean((f - d) ** 2)), np.max(np.abs(d))))
# tail vs quadrature
rho = lambda y: np.log(y / (2 * np.pi)) / (2 * np.pi)
for xx in [14.13, 500.0, 0.5 * Y0, Y0 - 50, Y0 - 1.0, Y0 - 0.3]:
    g = lambda y: 0.5 * rho(y) * (1 / (xx - y) + 1 / (xx + y))
    q1 = quad(g, Y0, Y0 + 10, limit=400, epsabs=1e-13)[0] + quad(g, Y0 + 10, np.inf, limit=400, epsabs=1e-13)[0]
    print('tail x=%.2f closed %.12f quad %.12f' % (xx, tail_force(np.array([xx]), Y0)[0], q1))
# total force in the middle: mean should be about -pi/8
tot = 0.5 * f + tail_force(x, Y0)
i = slice(int(0.2 * M), int(0.8 * M))
print('mean total force middle 60%%: %.5f (pi/8 = %.5f); 20 lowest / highest zeros mean force %.4f / %.4f'
      % (tot[i].mean(), np.pi / 8, tot[:20].mean(), tot[-20:].mean()))
