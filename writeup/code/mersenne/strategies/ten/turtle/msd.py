"""Diffusion of the prime turtle ("left at every prime") against random models: mean squared displacement over
windows of w numbers, for the numbers up to X.   usage: msd.py X"""
import sys, math, random
import numpy as np
X = int(sys.argv[1])
sv = np.ones(X + 1, bool); sv[:2] = False
for i in range(2, int(X ** .5) + 1):
    if sv[i]: sv[i * i::i] = False
DX = np.array([1, 0, -1, 0]); DY = np.array([0, 1, 0, -1])
def path(turn):
    """turn: boolean array over n = 1..X (index n). The turtle steps, then turns left if turn[n]."""
    t = turn[1:X + 1].astype(np.int64)
    head = np.concatenate(([0], np.cumsum(t)[:-1])) % 4      # heading during step n
    x = np.cumsum(DX[head]); y = np.cumsum(DY[head])
    return x, y
def msd(x, y, w):
    dx = x[w:] - x[:-w]; dy = y[w:] - y[:-w]
    return float(np.mean(dx * dx + dy * dy)) / w
rng = np.random.default_rng(7)
n = np.arange(X + 1, dtype=float); p1 = np.zeros(X + 1); p1[3:] = 1 / np.log(n[3:]); p2 = np.zeros(X + 1); p2[3::2] = 2 / np.log(n[3::2])
res = {}
xp, yp = path(sv)
W = [w for w in (10, 30, 100, 300, 1000, 3000) if w < X // 5]
res['primes'] = [msd(xp, yp, w) for w in W]
for name, pr in (('cramer', p1), ('cramer_odd', p2)):
    vals = []
    for s in range(20):
        tr = rng.random(X + 1) < pr
        x, y = path(tr)
        vals.append([msd(x, y, w) for w in W])
    v = np.array(vals); res[name] = (v.mean(0).tolist(), v.std(0).tolist())
print('X =', X, 'windows', W)
print('primes      ', ' '.join('%7.3f' % v for v in res['primes']))
for name in ('cramer', 'cramer_odd'):
    m, s = res[name]; print('%-11s ' % name, ' '.join('%7.3f' % v for v in m), ' (spread over 20 runs: ' + ' '.join('%.3f' % v for v in s) + ')')
