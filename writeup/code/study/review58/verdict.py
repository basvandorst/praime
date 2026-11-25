import json, math, sys, glob
import numpy as np
sys.path.insert(0, '.')
from r1_offsets import theta, jk_mean, load
def rows(files):
    out = []
    for fn in files:
        z = np.load(fn); info = json.loads(str(z['info'])); g = z['g']; gg = g[(g > 0) & (g <= 150)]
        out.append((info['key'], int(round(info['Q'])), info['eps'][0], info['c'], np.exp(2j * theta(gg, info['Q'], info['mu'])).sum().real + info['eps'][0] * info['c'] * 150 / (2 * math.pi)))
    return out
R = rows(sorted(glob.glob('rep/*.npz'))); x = np.array([r[4] for r in R]); keys = [r[0] for r in R]
m, s = jk_mean(x); print('all 45: %.3f +- %.3f  median %.2f  sd %.2f' % (m, s, np.median(x), x.std(ddof=1)))
sel = np.array([k != 'R_N1403_13' for k in keys]); m2, s2 = jk_mean(x[sel]); print('without curve 13 (44): %.3f +- %.3f' % (m2, s2))
first6 = np.array(['_%d.' % i in fn for fn in sorted(glob.glob('rep/*.npz')) for i in [0]])
idx = np.array([int(k.split('_')[-1]) for k in keys]); un = idx >= 6
m3, s3 = jk_mean(x[un]); print('the 39 not seen when the note was written (idx >= 6): %.3f +- %.3f' % (m3, s3))
for H, v in (('H_real', 4.97), ('H_fluct 1.41', 1.41), ('H_fluct 1.25', 1.25)):
    print('  distance to %-13s %+.2f = %.2f se' % (H, m - v, (m - v) / s))
print('  midpoint 4.97/1.41 = %.2f, 4.97/1.25 = %.2f' % ((4.97 + 1.41) / 2, (4.97 + 1.25) / 2))
d = 4.97 - m; sd = math.hypot(1.27, s); print('rounds 151+152 minus replication: %.2f +- %.2f (%.2f sigma)' % (d, sd, d / sd))
# residual vs conductor and vs c (trend inside the box?)
Ns = np.array([r[1] for r in R]); eps = np.array([r[2] for r in R]); cc = np.array([r[3] for r in R])
print('corr(res, N) %.3f  corr(res, eps*c) %.3f  (n=45, 1/sqrt(n)=%.3f)' % (np.corrcoef(x, Ns)[0, 1], np.corrcoef(x, eps * cc)[0, 1], 1 / math.sqrt(45)))
# rank test: replication vs old residuals
old = [r for r in load() if r['cls'] == 'g2']
xo = np.array([np.exp(2j * theta(r['g'][(r['g'] > 0) & (r['g'] <= 150)], r['Q'], r['mu'])).sum().real + r['eps'] * r['c'] * 150 / (2 * math.pi) for r in old])
from scipy.stats import mannwhitneyu, ttest_ind
print('old n=%d mean %.2f; Mann-Whitney p=%.3f; Welch t p=%.3f' % (len(xo), xo.mean(), mannwhitneyu(xo, x).pvalue, ttest_ind(xo, x, equal_var=False).pvalue))
print('old N range', min(int(round(r['Q'])) for r in old), max(int(round(r['Q'])) for r in old), ' rep N range', Ns.min(), Ns.max())
