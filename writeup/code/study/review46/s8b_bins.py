"""Gram-weighted cross statistic by conductor ratio and by parity; residual structure."""
import sys, math
import numpy as np
sys.path.insert(0, __import__('os').path.dirname(__file__))
from mylib import *
sig = sys.argv[1]
S = np.load(OUT + '/s8_gramx_%s.npz' % sig); S2 = np.load(OUT + '/s2_pairs.npz'); S1 = np.load(OUT + '/s1.npz')
I, J = S2['I'], S2['J']; Ds = np.array(S1['Ds']); y = S['x'].real; p = S['pred_H']; yi = S['x'].imag
d = np.abs(np.log(np.abs(Ds[J]) / np.abs(Ds[I])))
def fit(m, lab):
    Z = np.column_stack([np.ones(m.sum()), p[m]]); c = np.linalg.lstsq(Z, y[m], rcond=None)[0]
    reps = []
    for f in range(len(Ds)):
        k = m & ~((I == f) | (J == f)); Zk = np.column_stack([np.ones(k.sum()), p[k]]); reps.append(np.linalg.lstsq(Zk, y[k], rcond=None)[0])
    reps = np.array(reps); nj = len(reps); se = np.sqrt((nj - 1) / nj * ((reps - reps.mean(0)) ** 2).sum(0))
    print('  %-34s n %5d  slope %.4f +- %.4f  icpt %+.4f +- %.4f' % (lab, m.sum(), c[1], se[1], c[0], se[0]))
fit(d == 0, '|D1| = |D2|')
for lo, hi in ((1e-9, 0.01), (0.01, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.3)):
    fit((d >= lo) & (d < hi), '%.2f <= |ln(q2/q1)| < %.2f' % (lo, hi))
s1 = np.sign(Ds[I]); s2 = np.sign(Ds[J])
fit((s1 > 0) & (s2 > 0), 'both D > 0'); fit(s1 != s2, 'opposite signs'); fit((s1 < 0) & (s2 < 0), 'both D < 0')
ev = (Ds[I] % 2 == 0) | (Ds[J] % 2 == 0)
fit(ev, 'at least one D even'); fit(~ev, 'both D odd')
print('  Im X: mean %+.4f sd %.4f; corr(Im X, pred) %.3f' % (yi.mean(), yi.std(), np.corrcoef(yi, p)[0, 1]))
