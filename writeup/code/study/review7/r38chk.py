import json, math, numpy as np
from scipy.stats import kstest, gamma
d = json.load(open('/home/bas/praime/writeup/data/study/r38.json')); u = json.load(open('/home/bas/praime/writeup/data/study/r38_u.json'))
kl = lambda g: 'six' if g % 6 == 0 else ('ten' if g % 10 == 0 else 'other')
rat = {int(g): np.array(v['ratio5_4']) for g, v in d.items()}
for k in ('six', 'ten', 'other'):
    gs = [g for g in rat if kl(g) == k]
    f = np.polyfit(gs, [rat[g][0] for g in gs], 1)
    print(k, 'computed', gs, 'fit at lam16: g=250 %.3f g=300 %.3f' % tuple(np.polyval(f, [250, 300])))
g = np.array(u['g'])
for key in ('U3', 'U4', 'U5', 'U5_nearest', 'U_geometric'):
    v = np.array(u[key]); n = len(v)
    zs = (v.mean() - 1) / (v.std(ddof=1) / math.sqrt(n)); z0 = (v.mean() - 1) * math.sqrt(n)
    p_gamma = gamma.sf(v.sum(), n)
    print('%-11s mean %.3f se_sample %.3f z_sample %.2f | null se %.3f z_null %.2f one-sided p(Gamma) %.3f two-sided %.3f | KS %.3f'
          % (key, v.mean(), v.std(ddof=1) / math.sqrt(n), zs, 1 / math.sqrt(n), z0, p_gamma, 2 * p_gamma, kstest(v, 'expon').pvalue))
    for nm, m in (('six', g % 6 == 0), ('not6', g % 6 != 0)):
        w = v[m]; print('     %s n=%d mean %.3f se %.3f null-se %.3f' % (nm, len(w), w.mean(), w.std(ddof=1) / math.sqrt(len(w)), 1 / math.sqrt(len(w))))
six = g % 6 == 0
comp = np.array([s == 'computed' for s in u['order5']])
print('six computed', int((six & comp).sum()), 'of', int(six.sum()), '; others computed', int((~six & comp).sum()), 'of', int((~six).sum()))
