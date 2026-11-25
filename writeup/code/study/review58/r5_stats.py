"""Statistics of the genus-2 residual: window increments, robust location, regression models (additive vs c-scale),
and the same per unit T for Sym^3 and RS; noise per curve against the registered expectation."""
import json, math, sys
import numpy as np
sys.path.insert(0, '.')
from r1_offsets import theta, load, jk_mean
rows = load()
for r in rows:
    g = r['g']; g = g[(g > 0) & (g <= r['T2'])]
    r['gg'] = g; r['ph'] = np.exp(2j * theta(g, r['Q'], r['mu']))
def res(r, a, b):
    m = (r['gg'] > a) & (r['gg'] <= b)
    return r['ph'][m].sum() + r['eps'] * r['c'] * (b - a) / (2 * math.pi)
G2 = [r for r in rows if r['cls'] == 'g2']
print('genus 2 (55): window increments of Re(G - pred), jackknife; and t-test-like z')
for a, b in ((0, 10), (10, 30), (30, 60), (60, 90), (90, 120), (120, 150), (30, 150), (20, 150), (0, 150)):
    x = np.array([res(r, a, b).real for r in G2]); m, s = jk_mean(x)
    xp = np.array([res(r, a, b).real for r in G2 if r['eps'] > 0]); xm = np.array([res(r, a, b).real for r in G2 if r['eps'] < 0])
    print('  (%3d,%3d]  all %+5.2f +- %4.2f (z %+.1f)  per unit T/2pi %+.3f   eps+1 %+5.2f +- %4.2f  eps-1 %+5.2f +- %4.2f   sd/curve %.2f  median %+.2f' % (
        a, b, m, s, m / s, m / ((b - a) / 2 / math.pi), xp.mean(), jk_mean(xp)[1], xm.mean(), jk_mean(xm)[1], x.std(ddof=1), np.median(x)))
x = np.array([res(r, 0, 150).real for r in G2])
print('  sign test (0,150]: positive %d of %d;  trimmed mean (drop 5 each side) %+.2f' % ((x > 0).sum(), len(x), np.sort(x)[5:-5].mean()))
x2 = np.array([res(r, 30, 150).real for r in G2])
print('  sign test (30,150]: positive %d of %d' % ((x2 > 0).sum(), len(x2)))
# regression models at T = 150 (all genus 2 have T2 = 150): res = a + b * (eps c)
U = 150 / (2 * math.pi)
y = x; ec = np.array([r['eps'] * r['c'] for r in G2]); ep = np.array([r['eps'] for r in G2]); c = np.array([r['c'] for r in G2])
def ols(X, y):
    b = np.linalg.lstsq(X, y, rcond=None)[0]; n = len(y)
    loo = np.array([np.linalg.lstsq(np.delete(X, i, 0), np.delete(y, i), rcond=None)[0] for i in range(n)])
    se = np.sqrt((n - 1) / n * np.sum((loo - loo.mean(0)) ** 2, 0)); rss = np.sum((y - X @ b) ** 2)
    return b, se, rss
for name, X in (('additive a', np.c_[np.ones_like(y)]), ('c-scale b*(eps c U)', np.c_[ec * U]), ('a + b*eps c U', np.c_[np.ones_like(y), ec * U]),
                ('a + d*eps', np.c_[np.ones_like(y), ep]), ('a + g*c', np.c_[np.ones_like(y), c])):
    b, se, rss = ols(X, y)
    print('  model %-22s coef %s se %s  rss %.0f  (AIC-ish n ln(rss/n) + 2k = %.1f)' % (name, np.round(b, 3), np.round(se, 3), rss, len(y) * math.log(rss / len(y)) + 2 * X.shape[1]))
# per class, per unit T/2pi
for cls in ('g2', 'sym3', 'rs'):
    S = [r for r in rows if r['cls'] == cls]
    z = np.array([res(r, 0, r['T2']).real / (r['T2'] / 2 / math.pi) for r in S]); m, s = jk_mean(z)
    z2 = np.array([res(r, r['T2'] / 5, r['T2']).real / (0.8 * r['T2'] / 2 / math.pi) for r in S]); m2, s2 = jk_mean(z2)
    print('%-5s n=%2d  Re(G-pred)/(T/2pi) %+.3f +- %.3f ; on (T/5, T] %+.3f +- %.3f' % (cls, len(S), m, s, m2, s2))
# noise per curve against registration (Var/M along 0.152 genus 2; 0.479 Sym^3)
P = json.load(open('/home/bas/praime/writeup/code/study/lowzero/universal/blind4/predictions_blind4.json'))['L']
for cls in ('g2', 'sym3', 'rs'):
    S = [r for r in rows if r['cls'] == cls and r['rnd'] == 152]
    zz = np.array([res(r, 0, r['T2']) / P[r['key']]['sd']['along'] for r in S])
    print('%-5s round 152: Re(res)/sd_registered: mean %+.3f, rms %.3f (1 expected);  Im rms %.3f' % (cls, zz.real.mean(), math.sqrt(np.mean((zz.real - zz.real.mean()) ** 2)), math.sqrt(np.mean(zz.imag ** 2))))
