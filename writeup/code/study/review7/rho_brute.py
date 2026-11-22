"""Brute-force check of rho(k), k <= 3000: count nu_p directly from the residues of the point sets for every prime
11 <= p <= 7e5 (covers every prime factor of 210k +- g), compare with rho_mine (sieve over all primes).
The ratio brute/mine must be one constant (the generic tail beyond 7e5) for all k."""
import numpy as np, sys
sys.path.insert(0, '/home/bas/scratch/study-review7')
from rho_mine import primes_upto

B = 700_000
P = primes_upto(B); P = P[P >= 11].astype(np.int64); pf = P.astype(float)
K = np.arange(1, 3001)
rYY = np.load('/home/bas/scratch/study-review7/rhoYY_train.npy')


def nu(points):
    """points: list of int arrays (len K) -> nu_p for every (k, p): number of distinct residues."""
    R = np.stack([np.mod(x[:, None], P[None, :]) for x in points])   # (npts, K, nprimes)
    R.sort(axis=0)
    return 1 + (np.diff(R, axis=0) != 0).sum(0)


def S(points, chunk=40):
    out = np.zeros(len(points[0]))
    for i in range(0, len(points[0]), chunk):
        pts = [x[i:i + chunk] for x in points]
        v = nu(pts)
        out[i:i + chunk] = np.sum(np.log(1 - v / pf[None, :]) - len(points) * np.log(1 - 1 / pf[None, :]), axis=1)
    return out


for g in (2, 4, 6):
    d = 210 * K
    z = np.zeros_like(K); gg = np.full_like(K, g)
    Sg = S([z, gg])
    res = {'XX': S([z, gg, d, d + gg]) - 2 * Sg, 'XYp': S([z, gg, d]) - Sg, 'XYm': S([-d, z, gg]) - Sg}
    mine = np.load('/home/bas/scratch/study-review7/rho_g%d_train.npy' % g)
    for i, kind in enumerate(('XX', 'XYp', 'XYm')):
        ratio = np.exp(res[kind]) / mine[i][1:]
        print(g, kind, 'ratio brute/mine: min %.12f max %.12f spread %.2e' % (ratio.min(), ratio.max(), ratio.max() / ratio.min() - 1))
d = 210 * K; z = np.zeros_like(K)
ratio = np.exp(S([z, d])) / rYY[1:]
print('YY ratio min %.12f max %.12f spread %.2e' % (ratio.min(), ratio.max(), ratio.max() / ratio.min() - 1))
