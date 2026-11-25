import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+"/")
import math, time, numpy as np
from smlib import *
from scipy.special import jv
k = 4
rng = np.random.default_rng(0)
for (smin, smax) in [(1.8, 21.0), (6.0, 70.0)]:
    for mm in [1, 2, 5, 20, 100, 1000]:
        for fac in [0.35, 0.5]:
            nc = int(fac*(smax-smin)*smax/(math.pi*mm)) + 30
            nodes = cheb_nodes(smin, smax, nc)
            lam, E = kernel_factor(k, mm, nodes, tol=1e-13)
            s = rng.uniform(smin, smax, 300); t = rng.uniform(smin, smax, 300)
            Ls = bary_matrix(nodes, s); Lt = bary_matrix(nodes, t)
            approx = np.einsum('ir,r,ir->i', Ls @ E, lam, Lt @ E)
            ex = jv(k-1, s*t/(math.pi*mm))
            print('s in [%g,%g] m=%5d nc=%4d rank=%4d  max err %.1e (max|J| %.2e)' % (smin, smax, mm, nc, len(lam), np.max(np.abs(approx-ex)), np.max(np.abs(ex))))
