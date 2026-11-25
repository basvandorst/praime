"""S(t) = N(t) - theta(t)/pi on 600..1000 for every set of C and Bp (stored zeros), at the extremes (just before and
after each zero); mean of S over the band (Turing-type check: a missed zero shifts the mean by about -1 times the
fraction of the band above it); largest unfolded gap. Then Hardy-function sign changes in the largest gap of D = -395."""
import os, numpy as np, sys
from scipy.special import loggamma
R = os.path.dirname(os.path.abspath(__file__))+'/../../../../'
res = []
for f in ['r134_C.npz', 'r138_Bp.npz']:
    Z = np.load(R + 'writeup/data/riemann/lzeros/' + f)
    for k in Z.files:
        D = int(k[1:]); q = abs(D); a = 0 if D > 0 else 1
        z = np.sort(Z[k]); z = z[z > 0]
        zb = z[(z >= 600) & (z <= 1000)]
        th = lambda t: t / 2 * np.log(q / np.pi) + np.imag(loggamma((0.5 + a + 1j * t) / 2))
        n0 = np.searchsorted(z, zb - 1e-9)   # N just before each zero
        Sb = n0 - th(zb) / np.pi; Sa = Sb + 1
        # uniform grid mean of S
        tg = np.linspace(600, 1000, 40001)
        Sg = np.searchsorted(z, tg) - th(tg) / np.pi
        gaps = np.diff(zb) * np.log(q * zb[:-1] / (2 * np.pi)) / (2 * np.pi)
        res.append((D, max(np.abs(Sb).max(), np.abs(Sa).max()), Sg.mean(), np.abs(Sg).max(), gaps.max()))
res = np.array(res)
print('sets', len(res))
print('max|S| over extremes: quantiles 50/90/99/max: %.2f %.2f %.2f %.2f' % tuple(np.quantile(res[:, 1], [.5, .9, .99, 1])))
print('sets with max|S| > 1.5 at the extremes: %d' % (res[:, 1] > 1.5).sum(), res[res[:, 1] > 1.5][:, 0].astype(int))
print('mean S over band: range %.3f .. %.3f (sd %.3f)' % (res[:, 2].min(), res[:, 2].max(), res[:, 2].std()))
print('largest unfolded gap: sets > 3.0:', res[res[:, 4] > 3.0][:, [0, 4]])
for D in (197, -395):
    r = res[res[:, 0] == D][0]
    print('D %d: max|S| %.3f  mean S %+.4f  max gap %.3f' % (D, r[1], r[2], r[4]))
