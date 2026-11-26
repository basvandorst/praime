import numpy as np, sys
sys.path.insert(0, '/tmp/claude-1000/prover_boot')
from lp_ac import run, bfun, cellint
x0 = 1.5; U = np.log(30); n = 240
r, edges, ts, bs, base, M = run(x0, U, n, 50.0, 60.0, 0.02)
d = r.x[:-1]; m = r.x[-1]
print("m =", m)
# where is the constraint tight?
slack = bs - base - M @ d - m
tight = ts[slack < 1e-7]
print("tight t:", np.round(tight[:60], 2), "... count", len(tight))
# mass of Lambda in each cell relative to e^u: print cells with big weight, in x variable
w = (1 + d) * (np.exp(edges[1:]) - np.exp(edges[:-1]))     # Lambda mass of cell
x = np.exp(0.5*(edges[1:] + edges[:-1]))
order = np.argsort(-w)
print("cells with largest Lambda mass (x, mass, mass/width_x):")
for j in order[:25]:
    print("  x=%.3f  mass=%.3f  ratio=%.2f" % (x[j], w[j], 1 + d[j]))
print("fraction of cells at d=-1 (empty):", np.mean(d < -1 + 1e-9))
# check far t
tt = np.arange(60, 3000, 0.05)
Ff = -cellint(tt, 0, np.log(x0)) + np.array([cellint(tt, edges[j], edges[j+1]) for j in range(n)]).T @ d
bb = np.array([bfun(t) for t in tt[::20]])
print("min over t in [60,3000] of b - F (coarse):", np.min(bb - Ff[::20]))
