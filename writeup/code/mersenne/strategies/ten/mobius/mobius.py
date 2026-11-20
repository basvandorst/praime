"""Moebius walker: M(x) = sum_{n <= x} mu(n) (Mertens function) against random walkers that take a step +1 or -1 with
equal chance on the same square-free n (and stand still elsewhere).  Measures max |M(x)| / sqrt(x) over x in [200, X]
(Mertens's conjecture |M(x)| < sqrt(x) is false (Odlyzko and te Riele 1985) but holds far beyond any direct count),
and for R random walkers the share that ever exceeds sqrt(x) in [200, X] and their median max ratio.
Prediction written before the run: max |M(x)|/sqrt(x) stays below 0.6 up to 10^7, while more than half of the random
walkers exceed sqrt(x) somewhere in [200, 10^7].   usage: mobius.py X R seed out.txt"""
import sys, math
import numpy as np
X, R, seed, outp = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
out = open(outp, 'w')
def pr(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')
mu = np.ones(X + 1, np.int8); mu[0] = 0
isp = np.ones(X + 1, bool); isp[:2] = False
for p in range(2, X + 1):
    if not isp[p]: continue
    if p * p <= X: isp[p * p::p] = False
    mu[p::p] *= -1
    if p * p <= X: mu[p * p::p * p] = 0
M = np.cumsum(mu, dtype=np.int64)
xs = np.arange(X + 1, dtype=np.float64); sq = np.sqrt(xs)
lo = 200
r = np.abs(M[lo:]) / sq[lo:]; k = int(np.argmax(r))
pr('X = %d: square-free n <= X: %d (6/pi^2 X = %.0f); M(X) = %d' % (X, int(np.count_nonzero(mu)), 6 / math.pi ** 2 * X, M[X]))
pr('max |M(x)|/sqrt(x) on [200, X] = %.4f at x = %d (M = %d)' % (r[k], k + lo, M[k + lo]))
for e in range(3, int(math.log10(X)) + 1):
    a, b = max(lo, 10 ** (e - 1)), 10 ** e
    rr = np.abs(M[a:b + 1]) / sq[a:b + 1]; kk = int(np.argmax(rr))
    pr('  [%d, %d]: max %.4f at x = %d; M(10^%d) = %d' % (a, b, rr[kk], kk + a, e, M[b]))
sf = np.nonzero(mu)[0]
rng = np.random.default_rng(seed)
mx = []; over = 0
for t in range(R):
    steps = rng.integers(0, 2, size=sf.size, dtype=np.int8) * 2 - 1
    W = np.cumsum(steps, dtype=np.int32)
    # walker value at x = W at the last square-free n <= x; the max ratio over x is attained at square-free n or just before the next one; check at square-free n (ratio largest there between steps since sqrt grows)
    sel = sf >= lo
    ratio = np.abs(W[sel]) / np.sqrt(sf[sel].astype(np.float64))
    m = float(ratio.max()); mx.append(m); over += (m >= 1)
mx = np.array(mx)
pr('random walkers on the same square-free n, %d runs (seed %d): share with max |W|/sqrt(x) >= 1 on [200, X]: %.3f; median max ratio %.3f; smallest %.3f' % (R, seed, over / R, float(np.median(mx)), float(mx.min())))
pr('share of random walkers whose max ratio is below that of M: %.3f' % float(np.mean(mx < r[k])))
