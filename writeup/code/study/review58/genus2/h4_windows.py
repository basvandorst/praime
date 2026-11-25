"""H4/H3: where in t is the genus-2 deficit? Slopes on windows (a, b] with prediction -eps c (b - a)/2pi (registered
estimator, jackknife over curves), per sample; the offset of the lowest zeros; Sym^3 and RS of rounds 151/152 on the same
relative windows for comparison. Old c throughout (H1 shows c is right to < 1 %)."""
import sys, os, math; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from g2data import load, slope
rows = load(('g2', 'sym3', 'rs'))
def win(r, a, b):
    m = (r['gg'] > a) & (r['gg'] <= b); return complex(r['ph'][m].sum())
def jkmean(x):
    x = np.asarray(x, float); n = len(x); loo = (x.sum() - x) / (n - 1); return x.mean(), math.sqrt((n - 1) / n * np.sum((loo - loo.mean()) ** 2))
G2 = [r for r in rows if r['cls'] == 'g2']
samples = {'all100': G2, 'r151+152': [r for r in G2 if r['rnd'] != 'rep'], 'rep45': [r for r in G2 if r['rnd'] == 'rep'],
           'eps+1': [r for r in G2 if r['eps'] > 0], 'eps-1': [r for r in G2 if r['eps'] < 0]}
print('genus 2: slope on windows (a, b], pred = -eps c (b - a)/2pi')
W = [(0, 150), (0, 10), (10, 150), (20, 150), (0, 50), (50, 100), (100, 150), (0, 75), (75, 150), (37.5, 150)]
for nm, S in samples.items():
    out = []
    for a, b in W:
        bb, se = slope([win(r, a, b) for r in S], [-r['eps'] * r['c'] * (b - a) / (2 * math.pi) for r in S])
        out.append('(%g,%g] %.3f+-%.3f' % (a, b, bb, se))
    print('  %-9s n=%3d ' % (nm, len(S)) + '  '.join(out))
print('genus 2: mean residual Re(G - pred) on windows, all 100 / eps+1 / eps-1')
for a, b in [(0, 5), (0, 10), (10, 20), (0, 20), (20, 150), (10, 150)]:
    res = {nm: jkmean([(win(r, a, b) - (-r['eps'] * r['c'] * (b - a) / (2 * math.pi))).real for r in S]) for nm, S in samples.items()}
    print('  (%g,%g]  ' % (a, b) + '  '.join('%s %+.2f+-%.2f' % (nm, *res[nm]) for nm in ('all100', 'eps+1', 'eps-1', 'r151+152', 'rep45')))
# the bias of the registered slope from a common additive offset a: b_bias = sum Re(a conj p)/sum|p|^2
for a0 in (1.0, 1.58):
    P = np.array([-r['eps'] * r['c'] * 150 / (2 * math.pi) for r in G2]); print('  slope bias from a common offset %+.2f on (0,150]: %+.4f' % (a0, (a0 * P).sum() / (P ** 2).sum()))
# zeros in (0,10]: how many, and the mean phase of the first zero
n10 = [int(((r['gg'] > 0) & (r['gg'] <= 10)).sum()) for r in G2]
print('  zeros in (0,10]: mean %.2f, range %d-%d; Re e^{2i theta} of first zero: mean %+.3f' % (np.mean(n10), min(n10), max(n10), np.mean([r['ph'][0].real for r in G2])))
# other classes: same relative windows (fractions of T2)
print('Sym^3 and RS (rounds 151/152): slope on relative windows (fraction of T2)')
for cls in ('sym3', 'rs'):
    S = [r for r in rows if r['cls'] == cls]
    out = []
    for fa, fb in [(0, 1), (0, 1 / 15), (1 / 15, 1), (0, 1 / 2), (1 / 2, 1)]:
        bb, se = slope([win(r, fa * r['T2'], fb * r['T2']) for r in S], [-r['eps'] * r['c'] * (fb - fa) * r['T2'] / (2 * math.pi) for r in S])
        out.append('(%.3g,%.3g]T %.3f+-%.3f' % (fa, fb, bb, se))
    res = jkmean([(win(r, 0, r['T2'] / 15) - (-r['eps'] * r['c'] * r['T2'] / 15 / (2 * math.pi))).real for r in S])
    print('  %-5s n=%3d ' % (cls, len(S)) + '  '.join(out) + '   offset (0,T/15] %+.2f+-%.2f' % res)
