import json, glob, numpy as np
act = [0.18550989745840515, 0.28190597602467216, 0.27544222904150706]
DS = ('0.1', '0.03', '0.01')
r98 = json.load(open('/home/bas/praime/writeup/data/study/r98.json'))['sets']
sets = {'r98_poisson': [[v[d] for d in DS] for k, v in r98.items() if k.startswith('poisson')],
        'r98_gue_pinned': [[v[d] for d in DS] for k, v in r98.items() if k.startswith('gue')]}
for f in sorted(glob.glob('/home/bas/scratch/study-review32/var_*.json')):
    name = f.split('var_')[1].rsplit('_', 1)[0]
    sets.setdefault(name, []).extend(json.load(open(f)).values())
for name, rows in sets.items():
    a = np.array(rows)
    for i, d in enumerate(DS):
        x = a[:, i]
        print('%-15s n=%2d d=%-4s mean %.4f sd %.4f median %.4f z %+.2f frac_below_actual %.2f' % (name, len(x), d, x.mean(), x.std(ddof=1), np.median(x), (act[i]-x.mean())/x.std(ddof=1), np.mean(x < act[i])))
