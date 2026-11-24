"""Quadrature check of the arithmetic P-part of k[3]: the production grid (9,216 nodes) against a finer one (54,400
nodes), same prime cutoff (B1 over p <= 1000 in both, to keep the fine run short). usage: kf3_grid_check.py"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
import tplib
LS = [0.5 * i for i in range(1, 17)]; PS = tplib.PrimeSums(1000)
res = []
for kw in (dict(), dict(WS=0.125, NS=10, WT=0.25, NT=10)):
    q = tplib.Quad(**kw); d = tplib.T3_parts(0 * q.T, q.T, q.S, 9.4, (2,), PS=PS)
    v = q.kf3(d['P'] + d['Q']); res.append(np.array([v[int(round(l / q.WS)) - 1] for l in LS]))
    print('%-45s ' % kw + ' '.join('%+.7f' % x for x in res[-1]))
print('max |difference| = %.1e' % np.abs(res[0] - res[1]).max())
