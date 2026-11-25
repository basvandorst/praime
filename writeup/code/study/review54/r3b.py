import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+"/")
import json, numpy as np
from r1_sums import theta, own_data, D
pred = json.load(open(D + 'predictions.json'))['Lfunctions']
L = own_data()
Ts = np.arange(400, 601, 20.0)
print('T:', Ts)
for key in ('D12','F16','F18','F20','F22','F26'):
    v = L[key]; P = pred[key]
    g = np.load(D + 'zeros/zeros_%s.npz' % key)['g']; g = np.sort(g[(g > 0)])
    e = np.exp(2j * theta(g, v['Q'], v['mu'])); eps = v['eps']; c = P['c']; rot = np.conj(-eps)
    R = [ (e[g <= t].sum() + eps * c * t / (2*np.pi)) * rot for t in Ts]
    # Hann window residual on (0, T2)
    print(key, 'Im:', ' '.join('%+5.1f' % r.imag for r in R))
    print(key, 'Re:', ' '.join('%+5.1f' % r.real for r in R))
