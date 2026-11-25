import numpy as np, sys
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
from analysis import *
Dsel = [D for D in Ds if D != 1153]; a = arrays(Dsel); n = len(Dsel)
gs2 = np.array(P['gue_sigma2']); gkf = np.array(P['gue_kf3'])
kf_data = a['k3'] - 3*a['m2'] + 2*a['m1']
law_kf = a['kf3']
tot = law_kf.mean(0) - gkf
two = -3*(a['sigma2'].mean(0) - gs2)
three = a['full'].mean(0) - gue
print('law k[3] - GUE k[3]           ', np.round(tot, 3).tolist())
print('  of which -3(Sig2_law - Sig2_GUE)', np.round(two, 3).tolist())
print('  of which kappa3_law - kappa3_GUE', np.round(three, 4).tolist())
dd = kf_data - law_kf; se = dd.std(0, ddof=1)/np.sqrt(n)
print('data k[3] - law k[3] x1e4      ', np.round(dd.mean(0)*1e4, 1).tolist())
print('   se x1e4                     ', np.round(se*1e4, 1).tolist())
# alternative 'k[3] of GUE plus two-point terms': kappa3 = GUE kappa3 with the law's Sigma2 -> k[3]_alt = gue_kappa3 - 3 Sig2_law + 2l
alt = gue[None, :] - 3*a['sigma2'] + 2*ls[None, :]
da = kf_data - alt
print('data k[3] - alt(k3=GUE, Sig2 law) x1e4', np.round(da.mean(0)*1e4, 1).tolist(), ' z:', np.round(da.mean(0)/se, 1).tolist())
# data Sigma2 vs law
ds = a['m2'] - a['sigma2']; print('data Sigma2 - law x1e4', np.round(ds.mean(0)*1e4, 1).tolist(), 'se', np.round(ds.std(0, ddof=1)/np.sqrt(n)*1e4, 1).tolist())
