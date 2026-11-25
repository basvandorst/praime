"""H4: sharp against smooth windows. Hann sin^2(pi g/T) (registered R7; prediction half the sharp one), Tukey windows
(flat middle, cosine edges of width f*T), and sharp bands of 30 inside (0,150]; all 100 genus-2 curves, eps+1 and samples;
Sym^3 and RS with the same relative windows."""
import sys, os, math; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from g2data import load, slope
rows = load(('g2', 'sym3', 'rs'))
def tukey(x, f):
    w = np.ones_like(x); e = x < f; w[e] = np.sin(math.pi * x[e] / (2 * f)) ** 2; e = x > 1 - f; w[e] = np.sin(math.pi * (1 - x[e]) / (2 * f)) ** 2
    return w
def wsum(r, kind, par=None):
    x = r['gg'] / r['T2']
    if kind == 'hann': w = np.sin(math.pi * x) ** 2; mass = 0.5
    elif kind == 'tukey': w = tukey(x, par); mass = 1 - par
    elif kind == 'band': a, b = par; w = ((x > a) & (x <= b)).astype(float); mass = b - a
    return complex((w * r['ph']).sum()), mass
groups = {'g2 all100': [r for r in rows if r['cls'] == 'g2'], 'g2 eps+1': [r for r in rows if r['cls'] == 'g2' and r['eps'] > 0],
          'g2 r151+152': [r for r in rows if r['cls'] == 'g2' and r['rnd'] != 'rep'], 'g2 rep45': [r for r in rows if r['cls'] == 'g2' and r['rnd'] == 'rep'],
          'sym3': [r for r in rows if r['cls'] == 'sym3'], 'rs': [r for r in rows if r['cls'] == 'rs']}
specs = [('sharp', ('band', (0, 1))), ('hann', ('hann', None)), ('tukey .1', ('tukey', 0.1)), ('tukey .25', ('tukey', 0.25)), ('tukey .5', ('tukey', 0.5))] + \
        [('band %.1f-%.1f' % (a, a + 0.2), ('band', (a, a + 0.2))) for a in (0, 0.2, 0.4, 0.6, 0.8)]
for gname, S in groups.items():
    out = []
    for nm, (kind, par) in specs:
        Gs, P = [], []
        for r in S:
            g, mass = wsum(r, kind, par); Gs.append(g); P.append(-r['eps'] * r['c'] * mass * r['T2'] / (2 * math.pi))
        b, se = slope(Gs, P); out.append('%s %.3f+-%.3f' % (nm, b, se))
    print('%-12s n=%3d | ' % (gname, len(S)) + ' | '.join(out))
