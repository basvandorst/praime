import numpy as np, json, sys
from scipy import stats
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
from analysis import *
Dsel = [D for D in Ds if D != 1153]; a = arrays(Dsel); n = len(Dsel)
# (b) Hotelling calibration: null bootstrap of centred residual rows, and Gaussian with the residual covariance
d = a['k3'] - a['full']; dc = d - d.mean(0); C = np.cov(d.T)
rng = np.random.default_rng(5)
def T2of(x):
    m = x.mean(0); Cx = np.cov(x.T); return len(x)*m @ np.linalg.solve(Cx, m)
Tb = np.array([T2of(dc[rng.integers(0, n, n)]) for _ in range(4000)])
Tg = np.array([T2of(rng.multivariate_normal(np.zeros(16), C, n)) for _ in range(4000)])
Fref = lambda T: 1 - stats.f.cdf(T*(n - 16)/((n - 1)*16), 16, n - 16)
crit = stats.f.ppf(0.999, 16, n - 16)*(n - 1)*16/(n - 16)
print('Hotelling: F-reference 0.1%% point %.1f; null-bootstrap: mean T2 %.2f (F-ref mean %.2f), P(T2 > crit) = %.4f, 99.9%% quantile %.1f' % (
    crit, Tb.mean(), 16*(n - 1)/(n - 18), (Tb > crit).mean(), np.quantile(Tb, 0.999)))
print('           Gaussian sim: mean %.2f, P(T2 > crit) = %.4f; KS of F-ref p-values vs U(0,1) (bootstrap): p = %.3f' % (
    Tg.mean(), (Tg > crit).mean(), stats.kstest(Fref(Tb), 'uniform').pvalue))
print('           observed T2 %.2f -> bootstrap p %.3f' % (readings(a)['T2'], (Tb >= readings(a)['T2']).mean()))
# residual non-normality
sk = stats.skew(d, axis=0); ku = stats.kurtosis(d, axis=0)
print('skew of residuals per l:', np.round(sk, 2)); print('excess kurtosis:', np.round(ku, 2))
# single level offset (GLS) and its se
o = np.ones(16); Cm = C/n; mo = o @ np.linalg.solve(Cm, d.mean(0))/(o @ np.linalg.solve(Cm, o)); so = 1/np.sqrt(o @ np.linalg.solve(Cm, o))
print('GLS common level offset data - law: %+.2e +- %.2e' % (mo, so))
# (f) chi2 p-values for p = 7, 11, 13 with Hotelling-type reference (dof n-7)
r = readings(a); dof = n - 7
for p in (2, 3, 5, 7, 11, 13):
    c = r['chi2'][p]['full']; F_ = c*(dof - 16 + 1)/(dof*16); pv = 1 - stats.f.cdf(F_, 16, dof - 16 + 1)
    print('p=%2d chi2 vs law %.1f  p-value %.3f | vs GUE %.1f vs Q0 %.1f' % (p, c, pv, r['chi2'][p]['GUE'], r['chi2'][p]['Q0']))
    i = PR.index(p) + 1
    print('      coef x1e4 data', np.round(r['B'][i]*1e4, 0).astype(int).tolist())
    print('      coef x1e4 law ', np.round(r['pf'][i]*1e4, 0).astype(int).tolist(), ' se ', np.round(r['se'][i]*1e4, 0).astype(int).tolist())
# (d) include D = 1153
a2 = arrays(list(Ds)); r2 = readings(a2)
print('(d) with D=1153: T2 %.2f p %.3f gA %.4f gQ %.4f goff %.4f g2 %.5f max|m| %.2e; m(l) x1e4' % (r2['T2'], r2['pT'], r2['gA'], r2['gQ'], r2['goff'], r2['g2'], np.max(np.abs(r2['m']))), np.round(r2['m']*1e4, 2))
i1153 = Ds.index(1153)
print('    D=1153 kappa3 - law x1e4:', np.round((np.array(dat[1153]['k3']) - np.array(F['full'][i1153]))*1e4, 0).astype(int).tolist())
print('    resid sd x1e4          :', np.round(d.std(0)*1e4, 0).astype(int).tolist())
