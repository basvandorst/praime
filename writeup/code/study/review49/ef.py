# What the explicit formula alone delivers: registered readings with the law replaced by
#   EF_X = D_X + G_X : orthogonal lines of all prime powers n <= X with their weights W_p(k) (no zero used)
#                      + the GUE form factor beyond tau_X = ln X/L (conjecture, universal, no arithmetic)
#   D_X alone        : lines only (theorem-level part)
#   DH + GH          : Berry-type: lines up to the Heisenberg frequency (X = N t^2/4 pi^2) + GUE saturation beyond
import json, numpy as np, sys
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
import analyse as A
from scipy import stats

ids = sorted(i for i in A.LAW if i in A.S2 and i not in A.EXCL)
S = A.meas(ids); T = A.arr(ids, 'full', 'full'); n = len(ids)
G = np.array(A.GUE)
pick = [0, 1, 3, 7, 11, 15, 19]
print('curves', n, ' l =', A.LS[pick])
print('\nmean over curves (l = 0.5,1,2,4,6,8,10):')
print('  data - GUE              ', A.fmt((S - G).mean(0)[pick]))
print('  law  - GUE              ', A.fmt((T - G).mean(0)[pick]))
print('  data - law              ', A.fmt((S - T).mean(0)[pick]))
for X in ['10', '30', '100', '300', '1000', '10000', '100000', '1000000', '10000000', 'H2', 'H']:
    DX = A.arr(ids, 'full', 'D' + X); GX = A.arr(ids, 'full', 'G' + X)
    print('X=%-9s data - D_X %s | data - (D_X+G_X) %s | law - (D_X+G_X) %s' % (X, A.fmt((S - DX).mean(0)[pick], '%+.3f'),
          A.fmt((S - DX - GX).mean(0)[pick]), A.fmt((T - DX - GX).mean(0)[pick])))
# between-curve spread of law - (D_X + G_X) (the curve-specific part beyond X)
print('\nbetween-curve sd (l = 0.5,1,2,4,6,8,10):')
print('  data                    ', A.fmt(S.std(0)[pick], '%.4f'))
print('  data - law (noise)      ', A.fmt((S - T).std(0)[pick], '%.4f'))
for X in ['30', '100', '300', '1000', '10000', 'H']:
    DX = A.arr(ids, 'full', 'D' + X); GX = A.arr(ids, 'full', 'G' + X)
    print('  law - EF_X, X=%-8s  ' % X, A.fmt((T - DX - GX).std(0)[pick], '%.4f'), '  data - EF_X', A.fmt((S - DX - GX).std(0)[pick], '%.4f'))

# registered readings with EF_X as the template
print('\nregistered readings with the template replaced (bootstrap se):')
dk1p2 = A.arr(ids, 'full', 'dk1p2'); dw1p2 = A.arr(ids, 'full', 'dw1p2'); dk1 = A.arr(ids, 'full', 'dk1'); deg1 = A.arr(ids, 'full', 'deg1')
nosym = A.arr(ids, 'full', 'nosym')
LR = (np.array([A.P[i]['full'] for i in ids])[:, 11] - G[11]) > 0.02
for name, TT in [('law (own)', T)] + [('EF_%s' % X, A.arr(ids, 'full', 'D' + X) + A.arr(ids, 'full', 'G' + X)) for X in ('30', '100', '300', '1000', 'H')]:
    gA = A.gam(S, TT, np.zeros_like(TT), 500); gd = A.gam(S, TT, deg1, 500)
    gk = A.gam(S, TT, TT + dk1p2, 500); gw = A.gam(S, TT, TT + dw1p2, 500); gk1 = A.gam(S, TT, TT + dk1, 500)
    m = (S - TT).mean(0); se = (S - TT).std(0, ddof=1) / np.sqrt(n)
    Dd = S - TT; C = np.cov(Dd.T); T2 = n * m @ np.linalg.solve(C, m); F = T2 * (n - 20) / ((n - 1) * 20)
    lr = (S[LR] - G).mean(0); lrl = (TT[LR] - G).mean(0)
    print('  %-10s gA %.4f(%.4f) gdeg1 %.4f(%.4f) gk1p2 %.3f(%.3f) gw1p2 %.4f(%.4f) gk1 %.3f(%.3f) | level max|m|/se %.1f, max|m| %.4f, Hotelling p %.2g | T7 l=4,6: data %+.4f %+.4f template %+.4f %+.4f' % (
        name, *gA, *gd, *gk, *gw, *gk1, np.abs(m / se).max(), np.abs(m).max(), stats.f.sf(F, 20, n - 20), lr[7], lr[11], lrl[7], lrl[11]))
# Sym2 reading: amplitude of the Sym2 ratio in the off-diagonal term
print('\n  gamma_Sym2 (law vs Sym2 ratio set to 1): %.3f +- %.3f' % A.gam(S, T, nosym, 1000))
tm = T - nosym
print('  template law - nosym: between-curve sd per l', A.fmt(tm.std(0)[pick], '%.5f'), ' mean', A.fmt(tm.mean(0)[pick], '%+.5f'))
