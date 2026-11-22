"""Re-integrate round 54's expected counts from the stored log-shares: linear vs cubic interpolation, GL quadrature;
z values; post hoc split multiples of 6; a sub-window [1e10,1.2e10) check."""
import json, math
import numpy as np, mpmath
from scipy.interpolate import CubicSpline
from scipy.stats import poisson
r = json.load(open('/home/bas/praime/writeup/data/study/r54.json'))
cnt = {int(k): v for k, v in json.load(open('/home/bas/scratch/study-review13/c_1e10_2e10.json'))['counts'].items()}
L0, L1 = 1e10, 2e10; A, B = math.log(L0), math.log(L1)
GS = list(range(230, 421, 2)); LF = np.linspace(A, B, 9); L3 = np.array([A, B])
norm = r['gaps'] / float(mpmath.li(L1) - mpmath.li(L0))
x, w = np.polynomial.legendre.leggauss(40)
def law(g, a=A, b=B, kind='cubic'):
    u = 0.5 * (b - a) * x + 0.5 * (b + a); wu = 0.5 * (b - a) * w
    o2 = r['o2'][str(g)]
    ls = (CubicSpline(LF, o2)(u) if kind == 'cubic' else np.interp(u, LF, o2)) + np.interp(u, L3, r['o3'][str(g)])
    return float(np.sum(wu * np.exp(ls) * np.exp(u) / u))
lin = {g: law(g, kind='lin') * norm for g in GS}; cub = {g: law(g) * norm for g in GS}
rep = r['law_per_gap']
print('max rel diff cubic-vs-linear %.2e, mine-vs-repo %.2e' % (
    max(abs(cub[g] / lin[g] - 1) for g in GS), max(abs(lin[g] / rep[str(g)] - 1) for g in GS)))
for name, (lo, hi), sel in [('230-278', (230, 278), None), ('280-328', (280, 328), None), ('330-420', (330, 420), None),
                            ('pooled', (230, 420), None), ('pooled 6|g', (230, 420), 6), ('pooled 6!|g', (230, 420), -6)]:
    gg = [g for g in GS if lo <= g <= hi and (sel is None or (sel == 6 and g % 6 == 0) or (sel == -6 and g % 6))]
    c = sum(cnt.get(g, 0) for g in gg); lw = sum(cub[g] for g in gg)
    print('%-12s count %5d law %8.1f ratio %.4f +- %.4f z %+.2f z(0.84) %+.2f  P(X<=c|0.84law) %.1e' % (
        name, c, lw, c / lw, math.sqrt(c) / lw, (c - lw) / math.sqrt(lw), (c - .84 * lw) / math.sqrt(.84 * lw),
        poisson.sf(c - 1, .84 * lw)))
# exact relative position of the sub-window: [1e10,1.2e10) needs its own gap count; use the share at local lam
