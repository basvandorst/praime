"""Off-diagonal ('line') terms of the recipe in a sharp window: the density of e^{2i theta} carries, besides -eps c/2pi, the terms
(A_r/2pi) r^{it} with A_r = -eps c rho(r)/sqrt(ab) (the registered line amplitudes, r = 2, 1/2, 3, 1/3 from the constants files).
Their window integrals O_L(a,b) = sum_r (A_r/2pi) (r^{ib} - r^{ia})/(i ln r) are bounded and Q-independent in t, so if eps*c*rho(r)
is coherent over the family they give a COMMON t-structure.  Test: window means of the observed residual vs window means of O_L."""
import json, math, sys
import numpy as np
sys.path.insert(0, '.')
from r1_offsets import theta, load, jk_mean
B4 = '/home/bas/praime/writeup/code/study/lowzero/universal/blind4/constants4/'
B1 = '/home/bas/praime/writeup/code/study/lowzero/universal/blind/constants/'
R = {'2': 2.0, '1/2': 0.5, '3': 3.0, '1/3': 1 / 3}
rows = [r for r in load() if r['cls'] == 'g2']
for r in rows:
    d = json.load(open((B4 if r['rnd'] == 152 else B1) + r['key'] + '.json'))
    r['A'] = {k: complex(*v) for k, v in d['lines'].items()}
def O(r, a, b):
    return sum(r['A'][k] / (2 * math.pi) * (R[k] ** (1j * b) - R[k] ** (1j * a)) / (1j * math.log(R[k])) for k in R)
def res(r, a, b):
    g = r['g']; m = (g > a) & (g <= b)
    return np.exp(2j * theta(g[m], r['Q'], r['mu'])).sum() + r['eps'] * r['c'] * (b - a) / (2 * math.pi)
edges = list(range(0, 151, 10))
print('window   obs Re (jk)        O_L Re mean   | obs Im (jk)       O_L Im mean')
Wo = []; Wl = []
for a, b in zip(edges[:-1], edges[1:]):
    x = np.array([res(r, a, b) for r in rows]); y = np.array([O(r, a, b) for r in rows])
    Wo.append(x.mean()); Wl.append(y.mean())
    print('(%3d,%3d] %+5.2f +- %.2f     %+5.2f        | %+5.2f +- %.2f     %+5.2f' % (a, b, x.real.mean(), jk_mean(x.real)[1], y.real.mean(), x.imag.mean(), jk_mean(x.imag)[1], y.imag.mean()))
Wo = np.array(Wo[1:]); Wl = np.array(Wl[1:])
print('corr over windows (10,150]: Re %.3f  Im %.3f' % (np.corrcoef(Wo.real, Wl.real)[0, 1], np.corrcoef(Wo.imag, Wl.imag)[0, 1]))
for T in (10, 30, 60, 90, 120, 150):
    y = np.array([O(r, 0, T) for r in rows]); x = np.array([res(r, 0, T) for r in rows])
    print('cumulative T=%3d: O_L mean %+.2f %+.2fi (jk %.2f %.2f)   observed %+.2f %+.2fi' % (T, y.real.mean(), y.imag.mean(), jk_mean(y.real)[1], jk_mean(y.imag)[1], x.real.mean(), x.imag.mean()))
# per-curve regression: does subtracting O_L reduce the scatter / the window over-dispersion?
for lab, f in (('raw', lambda r, a, b: res(r, a, b)), ('minus O_L', lambda r, a, b: res(r, a, b) - O(r, a, b))):
    for part in ('real', 'imag'):
        m = []; s = []
        for a, b in zip(edges[1:-1], edges[2:]):
            x = np.array([f(r, a, b) for r in rows]); x = x.real if part == 'real' else x.imag
            m.append(x.mean()); s.append(x.std(ddof=1) / math.sqrt(len(x)))
        m = np.array(m); s = np.array(s)
        print('%-10s %s: chi2 of 14 window means vs 0 = %.1f' % (lab, part, np.sum((m / s) ** 2)))
print('mean eps*c*lambda-type coherence: mean of A_2, A_1/2, A_3, A_1/3 over curves:', {k: np.round(np.mean([r['A'][k] for r in rows]), 3) for k in R})
