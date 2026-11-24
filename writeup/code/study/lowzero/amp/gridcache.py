"""Build the Bogomolny-Keating grid of bklib.py (a copy of the repo's frozen bklib) once and cache its raw tables,
so that later scripts can rebuild the splines in a second. usage: gridcache.py  -> grid.npz"""
import os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'bk'))
import bklib
from scipy.interpolate import CubicSpline

CACHE = os.path.join(HERE, 'grid.npz')


class CachedGrid(bklib.Grid):
    def __init__(self):
        z = np.load(CACHE)
        self.eps = z['eps']; self.eps_max = float(self.eps[-1])
        self.raw = {k: z[k] for k in ('Z2', 'Z2m1', 'regd', 'Phi', 'Phim1', 'd2lnPhiDiag')}
        self.sp = {}
        for k, v in self.raw.items():
            self.sp[k] = (CubicSpline(self.eps, v.real), CubicSpline(self.eps, v.imag)) if np.iscomplexobj(v) else CubicSpline(self.eps, v)


def get():
    if not os.path.exists(CACHE):
        t0 = time.time(); G = bklib.Grid(eps_max=7.6, h=0.0025, pmax=10 ** 6)
        np.savez(CACHE, eps=G.eps, **G.raw); print('grid built and cached in %.0fs' % (time.time() - t0))
    return CachedGrid()


if __name__ == '__main__':
    G = get()
    print('eps grid', G.eps[0], G.eps[-1], len(G.eps))
