"""cramer_ctrl.py -- interpretive control: in a Cramer-type model (odd n is 'prime' independently with probability 2/ln n,
no zeros at all) the same regression gives coefficient 2*S_model = 4 on the model's own prime error: the factor 2 is the
response of a product of two densities to the fluctuation of one density, whatever causes that fluctuation.
Window 1e6..2e8 (cells of the same grid), band [2, 6000]."""
import numpy as np
import tzlib
from tz_counts import edges_for

rng = np.random.default_rng(2026)
lo, hi = 10 ** 4, 2 * 10 ** 8
ks, X = edges_for(lo, hi)
nc = len(X)
th = np.zeros(nc); pr = np.zeros(nc)
CH = 1 << 22
s = 1
prev_flag = False
while s <= hi:
    n = s + 2 * np.arange(CH + 1, dtype=np.int64)
    ln = np.log(n.astype(float))
    f = rng.random(CH + 1) < np.minimum(1.0, 2.0 / ln)
    f[n < 3] = False
    f[-1] = False          # the extra element is redrawn as the first of the next chunk; pairs across chunks are dropped (negligible)
    p = n[:-1][f[:-1]]; lp = ln[:-1][f[:-1]]
    keep = p <= hi
    j = np.searchsorted(X, p[keep])
    th += np.bincount(j, weights=lp[keep], minlength=nc)
    g = f[:-1] & f[1:]
    q = n[:-1][g]; lq = ln[:-1][g]
    keep = q <= hi
    pr += np.bincount(np.searchsorted(X, q[keep]), weights=(lq * np.log(q + 2.0))[keep], minlength=nc)
    s += 2 * CH
Xf = X.astype(float); u = np.log(Xf); sq = np.sqrt(Xf)
Yth = (np.cumsum(th) - Xf) / sq
# the model's own mean: sum over odd n<=x of 2/ln n * ln n = x (approximately); pairs: sum (2/ln n)(2/ln(n+2)) ln n ln(n+2) over odd n = 2x
Ypr = (np.cumsum(pr) - 2.0 * Xf) / sq
j0, j1 = int(np.searchsorted(X, 10 ** 6)), int(np.searchsorted(X, hi, side="right")) - 1
for glo, ghi in ((2, 300), (2, 6000)):
    b, cov, d = tzlib.fit(Ypr, [Yth], u, j0, j1, glo, ghi, 2)
    print("Cramer model, band [%g, %g]: coefficient of the pseudo-pair error on the pseudo-prime error %.3f +- %.3f (2*S_model = 4)" % (glo, ghi, b[0], np.sqrt(cov[0, 0])))
