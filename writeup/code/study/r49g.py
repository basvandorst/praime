"""Round 49 with direct quadrature (after the twelfth review): order 2 at 40 Gauss-Legendre nodes in u over
[ln 1e7, ln 1e8], the order-3 + overlap correction at three lam values with quadratic interpolation. usage: r49g.py"""
import os, json, math
import numpy as np
import mpmath
import r47
from r51 import o3_delta

r49 = json.load(open(os.path.join(r47.OUT, 'r49.json')))
L0, L1 = 10 ** 7, 10 ** 8
A, B = math.log(L0), math.log(L1)
XN, WN = np.polynomial.legendre.leggauss(40)
UN = 0.5 * (B - A) * XN + 0.5 * (B + A); WU = 0.5 * (B - A) * WN
L3 = np.array([A, 0.5 * (A + B), B])
norm = r49['gaps'] / float(mpmath.li(L1) - mpmath.li(L0))
law = {}
for g in range(92, 221, 2):
    st = r47.setup(g)
    ls2 = np.array([math.log(r47.share(st, float(u), 2)) for u in UN])
    d3 = np.polyval(np.polyfit(L3, o3_delta(g, L3), 2), UN)
    law[g] = float(np.sum(WU * np.exp(ls2 + d3) * np.exp(UN) / UN)) * norm
out = {}
for lo, hi in ((92, 118), (120, 148), (150, 178), (180, 220)):
    c = sum(r49['per_gap'][str(g)]['count'] for g in range(lo, hi + 1, 2)); lw = sum(law[g] for g in range(lo, hi + 1, 2))
    out['%d-%d' % (lo, hi)] = {'count': c, 'law': lw, 'ratio': c / lw, 'z': (c - lw) / math.sqrt(lw)}
    print(lo, hi, c, round(lw, 1), round(c / lw, 3), round((c - lw) / math.sqrt(lw), 2), flush=True)
json.dump({'bands': out, 'law_per_gap': {str(g): v for g, v in law.items()}}, open(os.path.join(r47.OUT, 'r49g.json'), 'w'), indent=1)
