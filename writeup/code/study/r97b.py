"""Round 97 diagnostic (after the run, not registered): the swing's mean square is (4/d^2) * 2Q (the factor 2: the pairs
(g_j, g_k) and (-g_j, -g_k) share the frequency and add coherently). The repulsion effect measured directly: Q for the
actual 10,000 zeros against Q for Poisson sets with the same mean count N(t) = (t/2pi) ln(t/2pi e) + 7/8 on the same
range (5 realisations, seeds 0-4), at d = 0.1. usage: r97b.py"""
import os, json, math
import numpy as np
from r97 import Z, q_actual

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
d = 0.1


def qset(g):
    G = np.concatenate([g, -g]); rho = 0.5 + 1j * G
    c2 = np.abs(((1 + d) ** rho - 1) / rho) ** 2
    import r97
    old = r97.G; r97.G = G
    try:
        return q_actual(c2)
    finally:
        r97.G = old


T = float(Z[-1])
t = np.linspace(2 * math.pi * math.e ** 0.0 + 1e-9, T, 2_000_000)
N = (t / (2 * math.pi)) * np.log(t / (2 * math.pi * math.e)) + 7 / 8
N = np.maximum.accumulate(N)
res = {'delta': d, 'Q_actual': qset(Z), 'Q_poisson': []}
for s in range(5):
    rng = np.random.default_rng(s)
    lo, hi = max(N[0], 0.0), N[-1]
    n = rng.poisson(hi - lo)
    pts = np.sort(np.interp(rng.uniform(lo, hi, n), N, t))
    res['Q_poisson'].append(qset(pts))
qa = res['Q_actual']; qp = np.array(res['Q_poisson'])
res['rms_actual'] = 2 / d * math.sqrt(2 * qa); res['rms_poisson_mean'] = float(np.mean(2 / d * np.sqrt(2 * qp)))
res['rms_poisson_sd'] = float(np.std(2 / d * np.sqrt(2 * qp)))
print(res)
json.dump(res, open(os.path.join(OUT, 'r97b.json'), 'w'), indent=1)
