"""Zero counting profile: mean over curves of n(t) - theta(t)/pi - c_exp (= S(t) + lower order) on a t grid, genus 2, Sym^3, RS;
per-curve contributions to the (30,60] window; and the residual profile in 10-unit windows."""
import json, math, sys
import numpy as np
sys.path.insert(0, '.')
from r1_offsets import theta, load, jk_mean
rows = load()
tg = np.arange(2, 150.01, 0.25)
for cls in ('g2', 'sym3', 'rs'):
    S = [r for r in rows if r['cls'] == cls]
    D = []
    for r in S:
        g = r['g']; g = g[(g > 0)]
        T2 = r['T2']; tt = tg[tg <= T2]
        n = np.searchsorted(g, tt, side='right')
        D.append(np.r_[n - theta(tt, r['Q'], r['mu']) / math.pi - r['info']['c_exp'], np.full(len(tg) - len(tt), np.nan)])
    D = np.array(D)
    m = np.nanmean(D, 0)
    print(cls, 'mean of N(t) - theta/pi - c_exp in t-bins of 10:')
    for a in range(0, 150, 10):
        sel = (tg > a) & (tg <= a + 10)
        v = D[:, sel]; pc = np.nanmean(v, 1)
        pc = pc[~np.isnan(pc)]
        if len(pc) < 3: continue
        mm, ss = jk_mean(pc)
        print('   (%3d,%3d]  %+.3f +- %.3f   (n=%d)' % (a, a + 10, mm, ss, len(pc)))
G2 = [r for r in rows if r['cls'] == 'g2']
print('\ngenus 2: residual increments in 10-unit windows, mean +- jk')
for a in range(0, 150, 10):
    x = []
    for r in G2:
        g = r['g']; m = (g > a) & (g <= a + 10)
        x.append((np.exp(2j * theta(g[m], r['Q'], r['mu'])).sum() + r['eps'] * r['c'] * 10 / (2 * math.pi)).real)
    mm, ss = jk_mean(x)
    print('   (%3d,%3d]  %+.2f +- %.2f' % (a, a + 10, mm, ss))
print('\nper-curve (30,60] residual, sorted')
v = []
for r in G2:
    g = r['g']; m = (g > 30) & (g <= 60)
    v.append(((np.exp(2j * theta(g[m], r['Q'], r['mu'])).sum() + r['eps'] * r['c'] * 30 / (2 * math.pi)).real, r['key'], r['eps']))
v.sort()
print(' '.join('%s:%+.1f' % (k, x) for x, k, e in v))
