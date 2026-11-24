"""RETROSPECTIVE. Split the per-L-function off-diagonal pattern x_i(l) = full_i - diag_i into the part carried by the
indicators of 2..13 (between groups: what the registered gamma sees) and the rest (within groups: other primes dividing
D, the non-additivity of the bad-prime product, the spread of L). Slope of the data on each part, with bootstrap errors,
in test B and in the rounds 117/120 sample (100-600).  usage: decomp.py -> decomp_out.txt"""
import os, sys, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'bk'))
from common import *

def parts(S, F, Dg, X, ix=None):
    ix = np.arange(len(S)) if ix is None else ix
    S, F, Dg, X = S[ix], F[ix], Dg[ix], X[ix]
    H = X @ np.linalg.pinv(X.T @ X) @ X.T
    x = F - Dg; y = S - Dg
    xb = H @ x; xw = x - xb                                  # between (fitted by the indicators) and within
    yb = H @ y; yw = y - yb
    xbc = xb - xb.mean(0); ybc = yb - yb.mean(0)
    gb = (xbc * ybc).sum() / (xbc ** 2).sum()
    gw = (xw * yw).sum() / (xw ** 2).sum()
    xc = x - x.mean(0); yc = y - y.mean(0)
    pl = (xc * yc).sum() / (xc ** 2).sum()
    share_w = (xw ** 2).sum() / (xc ** 2).sum()
    return gb, gw, pl, share_w

def run(nm, S, F, Dg, X):
    rng = np.random.default_rng(11); n = len(S)
    o = parts(S, F, Dg, X); bs = np.array([parts(S, F, Dg, X, rng.integers(0, n, n)) for _ in range(1000)])
    print('%-16s between (indicator) slope %.3f +- %.3f | within slope %.3f +- %.3f (share of x variance %.2f) | pooled %.3f +- %.3f | registered gamma %.3f' % (
        nm, o[0], bs[:, 0].std(), o[1], bs[:, 1].std(), o[3], o[2], bs[:, 2].std(), gamma(S, F, Dg, X)))
    return o, bs

Ds, S, F, Dg = load_B(); X = design(Ds)
run('B', S, F, Dg, X)
z = np.load(os.path.join(HERE, 'r120_100_600.npz')); D2 = list(z['D'])
run('r117/120 100-600', z['S'], z['F'], z['Dg'], design(D2))
for a, b in ((100, 300), (300, 600)):
    z = np.load(os.path.join(HERE, 'r120_%d_%d.npz' % (a, b)))
    run('r117/120 %d-%d' % (a, b), z['S'], z['F'], z['Dg'], design(D2))
# within B: what makes up the within part? the L-function groups by the set of primes <= 13 dividing D
x = F - Dg; H = X @ np.linalg.pinv(X.T @ X) @ X.T; xw = x - H @ x
big = [i for i, D in enumerate(Ds) if any(p > 13 for p in prime_factors(D))]
print('B: share of the within-part variance from sets with a prime > 13 dividing D: %.2f (%d of %d sets)' % ((xw[big] ** 2).sum() / (xw ** 2).sum(), len(big), len(Ds)))
