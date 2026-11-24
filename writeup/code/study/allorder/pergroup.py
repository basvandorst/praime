"""Per direction group, q = 5, k = 1, lam 22.69: slope per unit ln H (W(1e6) - W(1e4))/ln 100 of the profinite
CG(1e3) + CG(1e4) runs pooled (102 windows), against review44 dirs.json (registered windows, retrospective) and the
y = 17 linear response -V_d/lam of review44."""
import json, math, numpy as np
D = json.load(open('../../../../writeup/code/study/review44/dirs.json'))['5_1']
Do = json.load(open('../../../../writeup/code/study/review44/dirs_old.json'))['5_1']
lam = 22.6927
rows = {}
for fn in ('cgP_1000_22.69.json', 'cgP_10000_22.69.json'):
    for v in json.load(open(fn)).values():
        E = v['5']; lab = E['labels']; s = (np.array(E['blk']['1000000']) - np.array(E['blk']['10000'])) / math.log(100)
        for g in set(lab): rows.setdefault(g, []).append(s[[i for i, l in enumerate(lab) if l == g]].mean())
print('group            CG (102 windows)       registered (24 w)       old height (9 w)      -V_d/lam (y=17)')
for g in sorted(rows):
    x = np.array(rows[g]); m, e = x.mean(), x.std(ddof=1) / math.sqrt(len(x))
    print('%-15s  %+.4f +- %.4f      %+.4f +- %.4f       %+.4f +- %.4f      %+.4f' % (g, m, e, D[g]['slope'][0], D[g]['slope'][1], Do[g]['slope'][0], Do[g]['slope'][1], -D[g]['V'] / lam))
