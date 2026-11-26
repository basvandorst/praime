# Compute Xi(t) = xi(1/2+it) on a grid t in [0, TMAX], save to npz. mpmath, 30 digits.
import numpy as np, mpmath as mp, sys
mp.mp.dps = 30
def Xi(t):
    s = mp.mpf('0.5') + 1j*mp.mpf(t)
    v = 0.5*s*(s-1)*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
    return float(mp.re(v))
h = 0.01; TMAX = 70.0
ts = np.arange(0, TMAX+h/2, h)
xs = np.array([Xi(t) for t in ts])
np.savez('/tmp/claude-1000/conv_prover/xi_grid.npz', t=ts, xi=xs)
print(ts.size, xs[0], xs[-1], Xi(14.134725141734693))
