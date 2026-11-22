import numpy as np, math, sys, json, time
from gen import *
from q import Qparts
DS = (0.1, 0.03, 0.01)
def rms(g):
    _, same, opp, selfc = Qparts(np.sort(g))
    Q = 2*same + 2*opp
    return [2/d*math.sqrt(2*q) for d, q in zip(DS, Q)]
def p_g1(seed):
    rng = np.random.default_rng(1000+seed); n = rng.poisson(NT - N1)
    return Ninv(rng.uniform(N1, NT, n))
def p_palm(seed):
    return np.concatenate([[Z[0]], p_g1(seed)])
def g_stat(seed, lo=None):
    lo = N1 if lo is None else lo
    u, _ = gue_unfolded(2000+seed, radius='exact')
    rng = np.random.default_rng(3000+seed)
    o = u[10000] + rng.uniform(0, 1000)
    v = u - o; v = v[(v >= 0) & (v <= NT - lo)]
    return Ninv(lo + v)
def g_pin_exact(seed):
    u, _ = gue_unfolded(seed, radius='exact')
    v = u[10000:20000] - u[10000]
    return Ninv(N1 + v)
F = {'p_g1': p_g1, 'p_palm': p_palm, 'g_stat': g_stat, 'g_pin_exact': g_pin_exact, 'g_stat0': lambda s: g_stat(s, lo=0.0)}
mode, s0, s1 = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
out = {}
for s in range(s0, s1):
    t0 = time.time(); g = F[mode](s); out[s] = rms(g)
    print(mode, s, len(g), round(g[0], 3), [round(x, 4) for x in out[s]], round(time.time()-t0, 1), flush=True)
json.dump(out, open('/home/bas/scratch/study-review32/var_%s_%d.json' % (mode, s0), 'w'))
