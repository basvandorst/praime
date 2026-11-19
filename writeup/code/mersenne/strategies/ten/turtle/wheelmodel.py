"""Prediction made before the run: a Cramer model that respects the small primes (random numbers without the
factors in W, each with chance prod(p/(p-1)) / ln n) brings the turtle's diffusion and the gap spread close to the
primes. usage: wheelmodel.py X"""
import sys
import numpy as np
from math import gcd
X = int(sys.argv[1])
sv = np.ones(X + 1, bool); sv[:2] = False
for i in range(2, int(X ** .5) + 1):
    if sv[i]: sv[i * i::i] = False
DX = np.array([1, 0, -1, 0]); DY = np.array([0, 1, 0, -1])
def path(turn):
    t = turn[1:X + 1].astype(np.int64); head = np.concatenate(([0], np.cumsum(t)[:-1])) % 4
    return np.cumsum(DX[head]), np.cumsum(DY[head])
def msd(x, y, w):
    dx = x[w:] - x[:-w]; dy = y[w:] - y[:-w]; return float(np.mean(dx * dx + dy * dy)) / w
def gapstat(pr):
    P = np.nonzero(pr)[0]; P = P[P > 1000]; g = np.diff(P).astype(float); return (g * g).mean() / g.mean() - g.mean(), (g * g).mean() / g.mean() ** 2
rng = np.random.default_rng(11); n = np.arange(X + 1, dtype=float)
rows = [('primes', None)]
for W in ((), (2,), (2, 3), (2, 3, 5), (2, 3, 5, 7), (2, 3, 5, 7, 11, 13)):
    rows.append(('wheel ' + ('·'.join(map(str, W)) or 'none'), W))
print('X = %d; diffusion per step at windows 1000 and 3000; Var(g)/E[g]; E[g^2]/E[g]^2' % X)
for name, W in rows:
    if W is None:
        x, y = path(sv); vg, r = gapstat(sv); print('%-22s %7.3f %7.3f   %7.3f %6.3f' % (name, msd(x, y, 1000), msd(x, y, 3000), vg, r)); continue
    ok = np.ones(X + 1, bool)
    for p in W: ok[::p] = False
    c = 1.0
    for p in W: c *= p / (p - 1)
    pr = np.zeros(X + 1); pr[3:] = np.minimum(1, c / np.log(n[3:])); pr[~ok] = 0
    m1 = []; m3 = []; v = []; rr = []
    for s in range(8):
        tr = rng.random(X + 1) < pr; x, y = path(tr); m1.append(msd(x, y, 1000)); m3.append(msd(x, y, 3000)); a, b = gapstat(tr); v.append(a); rr.append(b)
    print('%-22s %7.3f %7.3f   %7.3f %6.3f   (spread of the diffusion over 8 runs %.3f)' % (name, np.mean(m1), np.mean(m3), np.mean(v), np.mean(rr), np.std(m1)))
