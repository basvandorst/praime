import json, math, os
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); R = json.load(open(D + '/s1_rows.json'))
cx = lambda v: complex(v[0], v[1])
var = {'T2': 0.192, 'T3': 0.330, 'T4a': 0.479, 'T4b': 0.152}
for c in ('T4a', 'T4b', 'T2'):
    for sgn in (1, -1):
        ks = [k for k in R if R[k]['cls'] == c and round(R[k]['eps'][0]) == sgn]
        p = np.array([cx(R[k]['pred']).real for k in ks]); m = np.array([cx(R[k]['G']).real for k in ks]); M = np.array([R[k]['M'] for k in ks])
        b = np.sum(m * p) / np.sum(p ** 2); se_meas = math.sqrt(np.sum(p ** 2 * var[c] * M)) / np.sum(p ** 2); se_mod = math.sqrt(np.sum(p ** 2 * 0.077 * M)) / np.sum(p ** 2)
        print('%-4s eps=%+d n=%2d slope %.3f  se model %.3f  se measured-noise %.3f  m0 %s' % (c, sgn, len(ks), b, se_mod, se_meas, sorted(set(R[k]['m0'] for k in ks))))
        if c == 'T4a' and sgn == -1:
            for k in ks: print('     %-24s c %.3f pred %+6.1f G %+6.1f%+6.1fi m0 %d' % (k, R[k]['c'], cx(R[k]['pred']).real, *R[k]['G'], R[k]['m0']))
