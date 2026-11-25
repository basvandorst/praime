"""Exploratory H3 variant: a truncated singular series. c_R = sum_{q <= R squarefree} mu(q) prod_{p | q} (1 - c_p),
c_p = G_p(1/p) (good p) or F_p(1/p)/(1 - 1/p) (bad p), from the own local data (cown, 51 curves). Does some R make
(c_R - 1) ~ 0.9 (c - 1), and does the slope of the zero sums on c_R reach 1?"""
import sys, os, math, glob; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from g2data import load, slope
rows = {r['key']: r for r in load()}
Rmax = 100000
spf = np.zeros(Rmax + 1, dtype=np.int64)
for i in range(2, Rmax + 1):
    if spf[i] == 0: spf[i::i][spf[i::i] == 0] = i
Rs = [3, 10, 30, 100, 300, 1000, 3000, 10000, 100000]
out = {}
for fn in sorted(glob.glob('cown/*.npz')):
    key = os.path.basename(fn)[:-4]; z = np.load(fn); p = z['p']; sg = z['sgn'] if 'sgn' in z.files else np.ones(len(p))
    cp = dict(zip(p.tolist(), (sg * np.exp(z['logG'])).tolist()))
    f = np.zeros(Rmax + 1); f[1] = 1.0
    for q in range(2, Rmax + 1):
        s = spf[q]; m = q // s
        f[q] = 0.0 if m % s == 0 else -f[m] * (1 - cp[int(s)])
    cs = np.cumsum(f); out[key] = [cs[R] for R in Rs]
keys = list(out)
print('relative (c_R - 1)/(c - 1) (median over the 37 eps=+1 curves with c > 1.2):')
pos = [k for k in keys if rows[k]['eps'] > 0 and rows[k]['c'] > 1.2]
for j, R in enumerate(Rs):
    rel = np.array([(out[k][j] - 1) / (rows[k]['c'] - 1) for k in pos]); print('  R = %6d: median %.3f, IQR %.3f..%.3f' % (R, np.median(rel), *np.percentile(rel, [25, 75])))
S = [rows[k] for k in keys]
for j, R in enumerate(Rs):
    b, se = slope([r['G'] for r in S], [-r['eps'] * out[r['key']][j] * 150 / (2 * math.pi) for r in S])
    print('  slope on c_R, R = %6d: %.3f +- %.3f' % (R, b, se))
b, se = slope([r['G'] for r in S], [-r['eps'] * r['c'] * 150 / (2 * math.pi) for r in S]); print('  slope on c (these 51): %.3f +- %.3f' % (b, se))
print('rss of Re G against b * pred with the slope fitted, 51 curves (does the per-curve shape prefer c_R over c?):')
y = np.array([r['G'].real for r in S]); U = np.array([-r['eps'] * 150 / (2 * math.pi) for r in S])
for j, R in enumerate(Rs + ['c']):
    x = U * (np.array([out[r['key']][j] for r in S]) if R != 'c' else np.array([r['c'] for r in S]))
    b = (x * y).sum() / (x * x).sum(); print('  %-7s rss %.0f (slope %.3f)' % (R, ((y - b * x) ** 2).sum(), b))
