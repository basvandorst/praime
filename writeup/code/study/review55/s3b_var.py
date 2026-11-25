import json, math, os
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); R = json.load(open(D + '/s1_rows.json'))
cx = lambda v: complex(v[0], v[1])
def est(keys, b):
    al = []; ac = []; alh = []
    for k in keys:
        r = R[k]; pr = cx(r['pred']); G = cx(r['G']); u = pr / abs(pr); M = r['M']
        d = G * u.conjugate() - b * abs(pr)
        al.append(d.real ** 2 / M); ac.append(d.imag ** 2 / M)
        prh = cx(r['pred_hann']); Gh = cx(r['Gh']); dh = Gh * u.conjugate() - b * abs(prh)
        alh.append(dh.real ** 2 / (3 * M / 8))
    return np.mean(al), np.mean(ac), np.mean(alh), len(keys)
reg = [k for k in R if R[k]['status'] == 'COMPLETE']
for c in ['T1', 'T2', 'T3', 'T4a', 'T4b']:
    ks = [k for k in reg if R[k]['cls'] == c]
    num = sum((cx(R[k]['G']) * cx(R[k]['pred']).conjugate()).real for k in ks); den = sum(abs(cx(R[k]['pred'])) ** 2 for k in ks)
    a, b_, h, n = est(ks, num / den)
    print('%-4s n=%2d  Var/M along %.3f (model 0.077)  across %.3f (model 0.21)  Hann along per 3M/8 %.3f' % (c, n, a, b_, h))
a, b_, h, n = est(reg, 1.0021); print('ALL  Var/M along %.3f across %.3f Hann %.3f' % (a, b_, h))
