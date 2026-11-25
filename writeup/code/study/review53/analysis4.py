"""gamma_no3 and slope with the inclusion-exclusion log control (z - sum ln r1); thirds with inclusion-exclusion."""
import math, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from stats import jk, slope, P, R, tri, pairs, TT, PP, SS, m
from analysis import gl, Pdens, MONOS3, MONOS2, MONOS1, EDGES
M1 = m['single_trap']; M2 = m['pair_trap']; M3 = m['triple_trap']
pf = np.array([TT[t]['full'] for t in tri])
rs = {D: M1[i] / SS[D]['P1'] for i, D in enumerate(R)}
rp = {p: M2[i] / PP[p]['full'] for i, p in enumerate(pairs)}
z = np.array([math.log(rp[(a, b)] * rp[(a, c)] * rp[(b, c)]) for a, b, c in tri])
zs = np.array([math.log(rs[a] * rs[b] * rs[c]) for a, b, c in tri])
lm = np.log(M3)
for c in ['no3', 'kirk', 'lead']:
    xc = np.log([TT[t][c] for t in tri]); xf = np.log(pf)
    g = jk(lambda k: slope(xf[k] - xc[k], lm[k] - z[k] + zs[k] - xc[k]), tri)
    print('gamma_%s incl-excl %.4f +- %.4f' % (c, *g))
print('slope incl-excl %.4f +- %.4f' % jk(lambda k: slope(np.log(pf[k]), lm[k] - z[k] + zs[k]), tri))
M3p = m['triple_parts']; M2p = m['pair_parts']; M1p = m['single_parts']
for ip, (a, b) in enumerate(EDGES[1:]):
    p3 = np.array([gl(lambda tt: Pdens(TT[t]['key'], TT[t]['hm'], tt, MONOS3), a, b) for t in tri])
    p2 = {p: gl(lambda tt: Pdens(p, PP[p]['hm'], tt, MONOS2), a, b) for p in pairs}
    p1 = {D: gl(lambda tt: Pdens((D,), SS[D]['hm1'], tt, MONOS1), a, b) for D in R}
    rpp = {p: M2p[i, ip] / p2[p] for i, p in enumerate(pairs)}; rss = {D: M1p[i, ip] / p1[D] for i, D in enumerate(R)}
    rr = M3p[:, ip] / p3
    ie = np.array([rr[i] - (rpp[(x, y)] + rpp[(x, w)] + rpp[(y, w)] - 3) + (rss[x] + rss[y] + rss[w] - 3) for i, (x, y, w) in enumerate(tri)])
    r2ie = np.array([rpp[p] - rss[p[0]] - rss[p[1]] + 2 for p in pairs])
    print('[%4.0f,%4.0f] triples incl-excl %.4f +- %.4f | pairs incl-excl %.4f +- %.4f' % (a, b, *jk(lambda k: ie[k].mean(), tri), *jk(lambda k: r2ie[k].mean(), pairs)))
