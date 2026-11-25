"""Levels of the alternatives under three estimators (raw, pair CV, inclusion-exclusion r3 - sum(r2-1) + sum(r1-1));
off-diagonal fraction lambda for triples and for pairs; three-body derivative fraction."""
import math, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from stats import jk, slope, P, R, tri, pairs, TT, PP, SS, m

M1 = m['single_trap']; M2 = m['pair_trap']; M3 = m['triple_trap']
pf = np.array([TT[t]['full'] for t in tri]); r3 = M3 / pf
rs = {D: M1[i] / SS[D]['P1'] for i, D in enumerate(R)}
rp = {p: M2[i] / PP[p]['full'] for i, p in enumerate(pairs)}
cs = np.array([rs[a] + rs[b] + rs[c] - 3 for a, b, c in tri])
cp = np.array([rp[(a, b)] + rp[(a, c)] + rp[(b, c)] - 3 for a, b, c in tri])
A3 = np.array([TT[t]['A3_over_Aij'] for t in tri])
alts = {'diag': np.array([TT[t]['diag'] for t in tri]), 'kirkA': np.array([TT[t]['kirk'] for t in tri]) * A3,
        'no3A': np.array([TT[t]['no3'] for t in tri]) * A3, 'no3': np.array([TT[t]['no3'] for t in tri]),
        'kirk': np.array([TT[t]['kirk'] for t in tri])}
est = {'raw': lambda k: r3[k], 'pairCV': lambda k: r3[k] - cp[k], 'incl-excl': lambda k: r3[k] - cp[k] + cs[k]}
for en, f in est.items():
    lvl = jk(lambda k: f(k).mean(), tri)
    line = '%-9s level vs full %.4f +- %.4f |' % (en, *lvl)
    for an, pc in alts.items():
        d = pc / pf
        frac = jk(lambda k: (d[k].mean() - f(k).mean()) / (d[k].mean() - 1), tri)
        line += ' %s: frac seen %.2f +- %.2f;' % (an, *frac)
    print(line)
# pairs: lambda for the pair diag (theorem case)
r2 = np.array([rp[p] for p in pairs]); d2 = np.array([PP[p]['diag'] / PP[p]['full'] for p in pairs])
print('pairs: mean P_diag/P_full %.4f ; fraction of (diag - full) seen %.2f +- %.2f' % (
    d2.mean(), *jk(lambda k: (d2[k].mean() - r2[k].mean()) / (d2[k].mean() - 1), pairs)))
r2s = np.array([rp[p] - rs[p[0]] - rs[p[1]] + 2 for p in pairs])
print('pairs incl-excl r2 - sum(r1 - 1): %.4f +- %.4f' % jk(lambda k: r2s[k].mean(), pairs))
