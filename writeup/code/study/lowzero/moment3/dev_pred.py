import os, sys, math, time, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import tmlib as T
DEV = [1001, -1003, 1005, -1048, 1068, -1195, 1253, -1267]
t0 = time.time()
tr = T.Triple((1001, -1003, 1005))
print('one triple %.1f s' % (time.time() - t0))
ts = np.array([100.0, 316.0, 1000.0])
print('h multilinear', dict(zip(tr.ml.monos, np.round(tr.hm, 6))))
print('g diag', tr.g)
print('P', tr.P(ts), 'lead', tr.lead(ts), 'diag', tr.diag(ts))
# convergence in PMAX
for pm in [10**4, 10**5, 10**6]:
    t1 = time.time(); tt = T.Triple((1001, -1003, 1005), pmax=pm)
    print(pm, tt.hm, tt.P(ts), '%.1fs' % (time.time() - t1))
