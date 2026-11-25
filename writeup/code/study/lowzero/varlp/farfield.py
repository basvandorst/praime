"""farfield.py: check of the first-order far-field law of the three-point change (derivation.md sect. 3.5):
for |v| >> 1 (one zero far from the pair 0, u),  dT3(0, u, v) -> -dT2_diag(v) * sin(2 pi u)/(pi u),
averaged over v in [v0, v0 + 1] (one period of the off-diagonal phase), against the ratios-conjecture T3 (P-part + Q-part, vlib = tplib formulas), p = 2 and 5, L = 10."""
import sys, os, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import vlib
from twod import dQ
L = 10.0; PS = vlib.PrimeSums(1000)
u = np.array([0.0, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]) + 1e-7
for p in (2, 5):
    for v0 in (25.3, 47.1, 80.6):
        xg, wg = np.polynomial.legendre.leggauss(16); dT = 0; model = 0
        for xx, ww in zip(xg, wg):                      # average over one period in v: removes the e^{2 pi i v} part
            v = np.full_like(u, v0 + (xx + 1) / 2); z = np.zeros_like(u)
            B1 = [PS.B1_all(a, b) for (a, b) in vlib.b1_args(z, u, v, L)]
            T = vlib.T3P(z, u, v, L, [(), (p,)], B1)
            dT = dT + ww / 2 * (T[1] - T[0] + dQ(L, p, z, u, v))
            model = model + ww / 2 * (-vlib.dT2_diag(L, p, v) * 2 * np.sinc(2 * u))
        print('p=%d v=%5.1f  u:      ' % (p, v0) + ' '.join('%8.2f' % x for x in u))
        print('    ratios dT3*1e3:   ' + ' '.join('%+8.4f' % x for x in 1e3 * dT))
        print('    model  -dC D*1e3: ' + ' '.join('%+8.4f' % x for x in 1e3 * model))
