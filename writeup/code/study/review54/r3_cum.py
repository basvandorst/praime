import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+"/")
import json, numpy as np
from r1_sums import theta, own_data, D
pred = json.load(open(D + 'predictions.json'))['Lfunctions']
L = own_data()
for key in ('X5','X5d','X7','X7d'):
    L[key] = dict(Q=float(pred[key]['Q']), mu=pred[key]['mu'], eps=complex(*pred[key]['eps']), fam='GL1')
for key, v in L.items():
    P = pred[key]; T2 = float(P['T2'])
    g = np.load(D + 'zeros/zeros_%s.npz' % key)['g']; g = np.sort(g[(g > 0) & (g <= T2)])
    e = np.exp(2j * theta(g, v['Q'], v['mu']))
    eps = v['eps']; c = P['c']
    # residual R(t) = (sum_{gamma<=t} e + eps c t/2pi) * conj(-eps)/|eps|  (rotated so that the main term is real positive)
    rot = np.conj(-eps) / abs(eps)
    ts = np.array([T2 * f for f in (0.05, 0.1, 0.25, 0.5, 0.75, 1.0)])
    R = []
    for t in ts:
        s = e[g <= t].sum() + eps * c * t / (2 * np.pi)
        R.append(s * rot)
    print('%-5s first zero %.2f | rotated residual at t/T2 = .05 .1 .25 .5 .75 1: %s' % (key, g[0], ' '.join('%+.1f%+.1fi' % (r.real, r.imag) for r in R)))
