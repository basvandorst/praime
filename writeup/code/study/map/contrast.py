"""mod-5 contrast: excess of the mixed (rc) directions minus the excess of the other five directions (eq, conj), per window"""
import json, math, numpy as np
def jk(v):
    v = np.asarray(v); n = len(v); loo = (v.sum() - v) / (n - 1); return math.sqrt((n - 1) / n * ((loo - loo.mean()) ** 2).sum())
out = {}
for name in ('new', 'old'):
    R = json.load(open('rep_%s.json' % name))
    for cell in ('5_1', '5_2', '7_1', '11_1'):
        T = R['cells'][cell]['types']
        if 'rc' not in T: continue
        others = [t for t in T if t != 'rc']; nd = {t: T[t]['ndir'] for t in T}
        s_rc = np.array(T['rc']['slope_pw']); s_ot = sum(nd[t] * np.array(T[t]['slope_pw']) for t in others) / sum(nd[t] for t in others)
        l_rc = T['rc']['lr']; l_ot = sum(nd[t] * T[t]['lr'] for t in others) / sum(nd[t] for t in others)
        c = (s_rc - l_rc) - (s_ot - l_ot)
        rr = s_rc.mean() / l_rc; ro = s_ot.mean() / l_ot
        print('%s q,k=%s: contrast (rc - rest) excess %+.4f +- %.4f (%+.1f);  slope/LR: rc %.2f, rest %.2f' % (name, cell, c.mean(), jk(c), c.mean() / jk(c), rr, ro))
        out['%s_%s' % (name, cell)] = [float(c.mean()), jk(c)]
a, b = out['new_5_1'], out['old_5_1']
w1, w2 = 1 / a[1] ** 2, 1 / b[1] ** 2; m = (w1 * a[0] + w2 * b[0]) / (w1 + w2); e = 1 / math.sqrt(w1 + w2)
print('mod 5 k=1 contrast combined: %+.4f +- %.4f (%+.1f)' % (m, e, m / e))
