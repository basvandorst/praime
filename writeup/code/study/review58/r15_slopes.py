"""Own registered-type slopes (sum Re(G conj pred)/sum |pred|^2, jackknife over L) per class / root number, round 152 alone and
with round 151's degree-4 sets; the Sym^3 watch-item rule; slopes after subtracting a constant offset per genus-2 curve."""
import math, sys
import numpy as np
sys.path.insert(0, '.')
from r1_offsets import theta, load
rows = load()
for r in rows:
    g = r['g']; g = g[(g > 0) & (g <= r['T2'])]
    r['G'] = np.exp(2j * theta(g, r['Q'], r['mu'])).sum(); r['pred'] = -r['eps'] * r['c'] * r['T2'] / (2 * math.pi)
def slope(S, off=0.0):
    X = np.array([((r['G'] - off).real * r['pred'], r['pred'] ** 2) for r in S]); sx, sy = X.sum(0); b = sx / sy
    loo = (sx - X[:, 0]) / (sy - X[:, 1]); n = len(S); se = math.sqrt((n - 1) / n * np.sum((loo - loo.mean()) ** 2)); return b, se, n
def show(lab, S, off=0.0):
    b, se, n = slope(S, off); print('%-44s n=%3d  b = %.3f +- %.3f' % (lab, n, b, se)); return b, se
for rnd in ((152,), (151,), (151, 152)):
    print('rounds', rnd)
    for cls in ('sym3', 'g2', 'rs'):
        for e in (None, 1, -1):
            S = [r for r in rows if r['rnd'] in rnd and r['cls'] == cls and (e is None or r['eps'] == e)]
            if len(S) >= 3: show('  %s eps %s' % (cls, e), S)
    S = [r for r in rows if r['rnd'] in rnd]; show('  all degree 4', S)
bm, sm, _ = slope([r for r in rows if r['rnd'] == 152 and r['cls'] == 'sym3' and r['eps'] < 0])
print('watch item (round 152 Sym3 eps=-1): |b-1| = %.3f <= 3 se = %.3f: %s ; b/se = %.1f >= 3: %s' % (abs(bm - 1), 3 * sm, abs(bm - 1) <= 3 * sm, bm / sm, bm / sm >= 3))
for off in (1.4, 5.0):
    print('genus 2 with a constant offset %.1f removed per curve:' % off)
    for rnd in ((152,), (151, 152)):
        for e in (None, 1, -1):
            S = [r for r in rows if r['rnd'] in rnd and r['cls'] == 'g2' and (e is None or r['eps'] == e)]
            show('  rounds %s eps %s' % (rnd, e), S, off)
    S = [r for r in rows if r['rnd'] == 152]
    X = np.array([((r['G'] - (off if r['cls'] == 'g2' else 0)).real * r['pred'], r['pred'] ** 2) for r in S]); sx, sy = X.sum(0)
    loo = (sx - X[:, 0]) / (sy - X[:, 1]); n = len(S)
    print('  all 120 (round 152) with genus-2 offset %.1f removed: b = %.3f +- %.3f' % (off, sx / sy, math.sqrt((n - 1) / n * np.sum((loo - loo.mean()) ** 2))))
