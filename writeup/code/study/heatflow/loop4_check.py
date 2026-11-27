"""Independent check of flowed zeros: H_t(2 gamma) = (1/8) E_G xi(1/2 - sqrt(t/2) G + i gamma), Gauss-Hermite + mpmath (no ODE)."""
import sys, numpy as np, mpmath as mp
mp.mp.dps = 25
nodes, weights = np.polynomial.hermite_e.hermegauss(60)   # weight exp(-g^2/2)
def xi(s):
    return mp.mpf(1)/2 * s * (s - 1) * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)
def Ht(gam, t):
    a = mp.sqrt(mp.mpf(t) / 2); tot = mp.mpf(0)
    for g, w in zip(nodes, weights):
        if g < 0: continue
        v = xi(mp.mpf(1)/2 + a * g + 1j * mp.mpf(gam))
        tot += w * (2 * v.real if g > 0 else v.real)
    return tot
d = np.load(sys.argv[1]); ts = d['t']; x = d['x']
for ti in [int(a) for a in sys.argv[2].split(',')]:
    t = ts[ti]; errs = []
    for k in [int(a) for a in sys.argv[3].split(',')]:
        g0 = x[ti, k]
        f = lambda g: Ht(g, t)
        try:
            r = mp.findroot(f, (mp.mpf(g0) - mp.mpf('0.02'), mp.mpf(g0) + mp.mpf('0.02')), solver='secant', tol=1e-12, maxsteps=30)
        except Exception as ex:
            print('  fail', k, ex); continue
        errs.append(float(r - g0))
        print(f"  t = {t:.2f}, zero {k}: flowed (ODE) {g0:.6f}, root of H_t {float(r):.6f}, diff {float(r-g0):+.2e}, moved from t = 0 by {g0 - x[0,k]:+.4f}")
    print(f"t = {t:.2f}: max |diff| {max(abs(e) for e in errs):.2e}")
