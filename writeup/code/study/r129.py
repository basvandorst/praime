"""Round 129 (preregister/round129.txt): zeta's number variance Sigma2(l) at L ~ 8.17 against Berry's sum over prime
powers (diagonal below tau = 1, saturation above) and against GUE (registered large-l form, and the exact sine kernel).
usage: r129.py"""
import os, json, math
import numpy as np
from scipy import integrate
from scipy.special import sici
from sympy import primerange
HERE = os.path.dirname(os.path.abspath(__file__)); RD = os.path.join(HERE, '..', '..', 'data', 'riemann')
g = np.array([float(l) for l in open(os.path.join(RD, 'zeros_40000_flint.txt')) if l.strip()])
Lz = np.log(g / (2 * math.pi)); m = (Lz >= 7.6) & (Lz <= 8.6)
u = (g / (2 * math.pi) * np.log(g / (2 * math.pi * math.e)) + 7 / 8)[m]; L = float(Lz[m].mean())
LS = [0.5 * i for i in range(1, 17)]
EG = 0.5772156649015329
def gue_reg(l): return (math.log(2 * math.pi * l) + EG + 1 - sici(2 * math.pi * l)[1]) / math.pi ** 2
def gue_exact(l): return l - 2 * integrate.quad(lambda x: (l - x) * (math.sin(math.pi * x) / (math.pi * x)) ** 2 if x else l, 0, l, limit=200)[0]
pp = [(p, k) for p in primerange(2, int(math.exp(L)) + 1) for k in range(1, 40) if p ** k <= math.exp(L)]
def berry(l):
    s = sum(math.sin(math.pi * l * k * math.log(p) / L) ** 2 / (k * k * p ** k) for p, k in pp)
    tail = integrate.quad(lambda t: math.sin(math.pi * l * t) ** 2 / t ** 2, 1, 200, limit=2000)[0] + 0.5 / 200
    return 2 / math.pi ** 2 * (s + tail)
rows = []
for l in LS:
    x = np.arange(u[0], u[-1] - l, 0.01); cnt = np.searchsorted(u, x + l) - np.searchsorted(u, x)
    rows.append({'l': l, 'data': float(np.var(cnt)), 'berry': berry(l), 'gue_registered': gue_reg(l), 'gue_exact': gue_exact(l)})
    r = rows[-1]; print('l %.1f data %.4f  Berry %.4f  GUE(reg) %.4f  GUE(exact) %.4f' % (l, r['data'], r['berry'], r['gue_registered'], r['gue_exact']))
r1 = all(abs(r['berry'] - r['data']) < abs(r['gue_registered'] - r['data']) for r in rows if r['l'] >= 2)
r2 = max(abs(r['berry'] - r['data']) for r in rows if r['l'] >= 1) < 0.02
r8 = [r for r in rows if r['l'] == 8.0][0]
print('mean L %.3f, n %d; readings %s %s; at l = 8: data - GUE(reg) %+.4f, data - GUE(exact) %+.4f' % (L, len(u), r1, r2, r8['data'] - r8['gue_registered'], r8['data'] - r8['gue_exact']))
json.dump({'L': L, 'n': int(len(u)), 'rows': rows, 'readings': [bool(r1), bool(r2)]}, open(os.path.join(HERE, '..', '..', 'data', 'study', 'r129.json'), 'w'), indent=1)
