import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+"/")
import json, math, numpy as np
from scipy.special import erf
from r1_sums import theta, own_data, D
pred = json.load(open(D + 'predictions.json'))['Lfunctions']
L = own_data()
for key in ('X5', 'X5d', 'X7', 'X7d'):
    L[key] = dict(Q=float(pred[key]['Q']), mu=pred[key]['mu'], eps=complex(*pred[key]['eps']), fam='GL1 complex character')
def F_gue(s): return erf(2 * s / math.sqrt(math.pi)) - (4 * s / math.pi) * np.exp(-4 * s * s / math.pi)
allg = []; fam = {}
for key, v in L.items():
    T2 = float(pred[key]['T2'])
    g = np.load(D + 'zeros/zeros_%s.npz' % key)['g']; g = np.sort(g[(g > 0) & (g <= T2)])
    u = theta(g, v['Q'], v['mu']) / math.pi       # unfolded
    s = np.diff(u)
    allg.append(s); fam.setdefault(v['fam'], []).append(s)
for name, ss in list(fam.items()) + [('ALL 20', allg)]:
    s = np.concatenate(ss); n = len(s)
    row = []
    for x in (0.1, 0.2, 0.3, 0.5):
        e = n * F_gue(x); o = (s < x).sum()
        row.append('<%.1f: %d obs / %.1f GUE' % (x, o, e))
    print('%-24s gaps %5d  mean %.4f  min %.4f | %s' % (name, n, s.mean(), s.min(), ' ; '.join(row)))
