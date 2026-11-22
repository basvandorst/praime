"""First-occurrence statistics from r24.json / r25.json and the raw OeS table: KS, median, dependence, trends."""
import json, math
import numpy as np
from scipy.stats import kstest, spearmanr, expon, pearsonr

r24 = json.load(open('/home/bas/praime/writeup/data/study/r24.json'))
r25 = json.load(open('/home/bas/praime/writeup/data/study/r25.json'))
g = np.array(r24['g']); P = np.array(r24['P'], float); U = np.array(r24['U']); U5 = np.array(r25['U'])
# raw table P(g)
Pt = {}
for line in open('/tmp/claude-1000/-home-bas-praime/cfd838fc-8861-4900-931a-268c5006ce39/scratchpad/dl_oes/t0.txt'):
    if line.startswith('#') or not line.strip() or line.startswith('-') or line.strip().startswith('g'):
        continue
    f = line.split()
    if '?' not in f[1]: Pt[int(f[0].rstrip('*'))] = int(f[1].rstrip('*'))
assert all(Pt[int(a)] == int(b) for a, b in zip(g, P)), 'P mismatch'
print('n', len(g), 'g range', g.min(), g.max(), 'P range %.3g %.3g' % (P.min(), P.max()), '456 in table:', 456 in Pt, 'P(456)=%s' % Pt.get(456))
lam = np.log(P); ratio = g / lam
print('g/ln P(g): min %.1f max %.1f' % (ratio.min(), ratio.max()))
for name, v in (('r24', U), ('r25', U5)):
    ks = kstest(v, 'expon')
    print(name, 'mean %.3f se %.3f median %.3f KS D %.4f p %.4f' % (v.mean(), v.std(ddof=1) / math.sqrt(len(v)), np.median(v), ks.statistic, ks.pvalue))
print('max |U25/U24 - 1| = %.3f; mean U ratio %.4f' % (np.abs(U5 / U - 1).max(), U5.mean() / U.mean()))
# a single overall scale c: U/c vs Exp(1)
for c in (1.0, 1.2, 1.27, 1.4):
    ks = kstest(U / c, 'expon'); print('scale c=%.2f: KS p %.3f median %.3f' % (c, ks.pvalue, np.median(U / c)))
# median sd for Exp(1) with n=138
print('sd of sample median of 138 Exp(1): %.3f' % (1 / (2 * 0.5 * math.sqrt(len(U)))))
# dependence: lag-1 correlations of U along g, of ln P
print('Spearman U_g vs U_g+2: r=%.3f p=%.3f' % spearmanr(U[:-1], U[1:]))
print('Spearman of expon-uniform F(U) lag1: %.3f' % spearmanr(1 - np.exp(-U[:-1]), 1 - np.exp(-U[1:]))[0])
# proximity of first occurrences: same region?
o = np.argsort(P); Ps = P[o]; gs = g[o]
d = np.diff(Ps)
print('closest pairs of first occurrences (distance, g1, g2):', sorted(zip(d, gs[:-1], gs[1:]))[:5])
print('pairs within 1e6: %d, within 1e7: %d' % ((d < 1e6).sum(), (d < 1e7).sum()))
# trend of U with lambda and with g/lambda
print('Spearman U vs ln P: r=%.3f p=%.3g;  U vs g/lnP: r=%.3f p=%.3g' % (*spearmanr(U, lam), *spearmanr(U, ratio)))
for lo, hi in ((16, 19), (19, 21), (21, 23), (23, 25), (25, 27.5)):
    m = (lam >= lo) & (lam < hi)
    v = U[m]
    print('ln P in [%g,%g): n %3d  mean U %.2f  median %.2f  mean ln U %.2f (Exp(1): -0.577, se %.2f)' % (
        lo, hi, m.sum(), v.mean(), np.median(v), np.log(v).mean(), 1.28 / math.sqrt(m.sum())))
for lo, hi in ((9, 12), (12, 14), (14, 16), (16, 19)):
    m = (ratio >= lo) & (ratio < hi); v = U[m]
    print('g/lnP in [%g,%g): n %3d  mean U %.2f  median %.2f  mean ln U %.2f' % (lo, hi, m.sum(), v.mean(), np.median(v), np.log(v).mean()))
# simulation under independent Exp(1) with rate scale c, KS p distribution: fraction of KS p<0.01 for c = 1
rng = np.random.default_rng(0)
print('P(KS p<0.01 | Exp(1), n=138) sim: %.4f' % np.mean([kstest(rng.exponential(size=138), 'expon').pvalue < 0.01 for _ in range(2000)]))
# record (maximal) gaps among the 138: P(g') > P(g) for all g' > g
rec = [int(x) for x in g if all(Pt[k] > Pt[int(x)] for k in Pt if k > x and k % 2 == 0)]
print('record gaps among them:', len(rec))
mrec = np.isin(g, rec)
print('U at record gaps: n %d mean %.2f median %.2f; non-record: mean %.2f median %.2f' % (mrec.sum(), U[mrec].mean(), np.median(U[mrec]), U[~mrec].mean(), np.median(U[~mrec])))
# bins by g (exogenous; bins by realized P(g) are confounded: an early first occurrence has small U AND small ln P)
from scipy.stats import chi2
print('--- by g (exogenous) ---')
for lo, hi in ((180, 250), (250, 320), (320, 390), (390, 456)):
    m = (g >= lo) & (g < hi); v = U[m]; k = m.sum()
    c = v.mean(); lo95, hi95 = 2 * k * c / chi2.ppf(0.975, 2 * k), 2 * k * c / chi2.ppf(0.025, 2 * k)
    print('g in [%d,%d): n %d  MLE scale c = mean U %.2f (95%% CI %.2f-%.2f)  median %.2f  typical ln P %.1f  g/lnP~%.1f' % (
        lo, hi, k, c, lo95, hi95, np.median(v), np.median(lam[m]), np.median(ratio[m])))
k = len(U); c = U.mean()
print('all: c = %.3f, 95%% CI %.2f-%.2f' % (c, 2 * k * c / chi2.ppf(0.975, 2 * k), 2 * k * c / chi2.ppf(0.025, 2 * k)))
print('Spearman U vs g: r=%.3f p=%.3f' % spearmanr(U, g))
