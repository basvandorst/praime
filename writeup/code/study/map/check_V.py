"""compare the review's V per group (dirs.json, dirs_old.json; responses pooled over the 24 windows) with map.py's
'rev' estimator (same counts, same noise correction) and with map.py's primary estimator"""
import json, sys, math, numpy as np
sys.path.insert(0, '.')
import map as M
R44 = '../../../../writeup/code/study/review44/'
LR = json.load(open('lr_rev24.json'))
for setname, ref in (('new', 'dirs.json'), ('old', 'dirs_old.json')):
    B = json.load(open(R44 + ref)); print('==', setname, ref)
    for q in (3, 5, 7):
        G = M.group(q); lr = LR['%s|%d_1' % (setname, q)]
        S = M.MSet('%s=val_%s' % (setname, setname)); _, dg, _, _ = S.basis(q, 1)
        Vrev = np.array(lr['rev']['V']) / lr['scale']          # the review stored V unscaled
        Vp = np.array(lr['y17']['V']) / lr['scale']
        mx = 0; mxp = 0
        for gi, g in enumerate(G.groups):
            v = B['%d_1' % q][g['label']]['V']; sel = dg == gi
            mx = max(mx, abs(Vrev[sel].mean() - v)); mxp = max(mxp, abs(Vp[sel].mean() / v - 1))
        print('  q=%d: %d groups, max |V_rev - V_review| = %.2e;  primary estimator (lag-0 cross products) differs by at most %.2f %%' % (q, len(G.groups), mx, 100 * mxp))
