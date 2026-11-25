import json, math, os
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); R = json.load(open(D + '/s1_rows.json'))
cx = lambda v: complex(v[0], v[1])
reg = [k for k in R if R[k]['status'] == 'COMPLETE']
def pc(keys):
    num = np.array([cx(R[k]['G']) * cx(R[k]['pred']).conjugate() for k in keys]); den = np.array([abs(cx(R[k]['pred'])) ** 2 for k in keys])
    s = num.sum() / den.sum(); n = len(keys)
    loo = np.array([(num.sum() - num[i]) / (den.sum() - den[i]) for i in range(n)])
    return s, math.sqrt((n - 1) / n * np.sum((loo.real - loo.real.mean()) ** 2)), math.sqrt((n - 1) / n * np.sum((loo.imag - loo.imag.mean()) ** 2))
neg = [k for k in reg if R[k]['eps'][0] < -0.5 and abs(R[k]['eps'][1]) < 1e-9]
print('D6 eps=-1 n=%d: slope %.4f jackknife se %.4f; classes %s' % ((len(neg),) + pc(neg)[:2] + (sorted(set(R[k]['cls'] for k in neg)),)))
for c in ('T2', 'T4a', 'T4b'):
    ks = [k for k in neg if R[k]['cls'] == c]; s = pc(ks); print('   %s n=%d slope %.3f ± %.3f' % (c, len(ks), s[0].real, s[1]))
cpx = [k for k in reg if abs(R[k]['eps'][1]) > 1e-9]
s = pc(cpx); print('D7 complex eps n=%d: s = %.4f%+.4fi, jackknife se re %.4f im %.4f' % (len(cpx), s[0].real, s[0].imag, s[1], s[2]))
allc = cpx + ['T1_17_5', 'T1_17_7']
for k in ('T1_17_5', 'T1_17_7'):
    z = json.load(open(D + '/s2_%s.json' % k)); R[k]['G'] = z['G']
s = pc(allc); print('D7 with chi mod 17 n=%d: s = %.4f%+.4fi, jackknife se re %.4f im %.4f' % (len(allc), s[0].real, s[0].imag, s[1], s[2]))
