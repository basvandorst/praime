"""oned2.py: the one-dimensional pieces (oned.py) on the L grid for every bad set of twod2.py (single primes and the
non-additive part of pairs). Off-diagonal integrals damped with E0 = 15 and 30 (raw); the two agree to < 1e-5 (oned_out).
usage: oned2.py OUT.json
"""
import sys, os, math, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import vlib
from oned import rnodes, closed
from twod2 import SINGLES, PAIRS
TWO_PI = vlib.TWO_PI
LS = (7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0)


def ints(h, r, w, e, E0):
    d = np.exp(-(e / E0) ** 2); lr = np.log(r)
    return dict(mu=2 * np.sum(w * lr * h * d), A=2 * np.sum(w * lr ** 2 * h * d), B=-(math.pi ** 2 / 2) * 2 * np.sum(w * r * h * d))


out = {'L': LS, 'singles': {}, 'pairs': {}}
for L in LS:
    Rmax = 150.0 * L / TWO_PI; r, w = rnodes(Rmax); e = TWO_PI * r / L
    for p in SINGLES:
        c = closed(p, L); h = vlib.dT2_off(L, p, r)
        a, b = ints(h, r, w, e, 15.0), ints(h, r, w, e, 30.0)
        out['singles'].setdefault(str(p), []).append(dict(L=L, **c, dmu_off=b['mu'], dA_off=b['A'], dB_off=b['B'],
                                                          conv=max(abs(a[k] - b[k]) for k in a)))
    for (p, q) in PAIRS:
        h = vlib.dT2_off(L, p, r, rest=(q,)) - vlib.dT2_off(L, p, r)
        b = ints(h, r, w, e, 30.0)
        out['pairs'].setdefault('%d,%d' % (p, q), []).append(dict(L=L, na_mu=b['mu'], na_A=b['A'], na_B=b['B']))
    print('L=%.1f done; p=2: dmu_off %+.5f dA_off %+.5f dB_off %+.5f; pair (2,3) NA mu %+.6f A %+.6f B %+.6f' % (
        L, out['singles']['2'][-1]['dmu_off'], out['singles']['2'][-1]['dA_off'], out['singles']['2'][-1]['dB_off'],
        out['pairs']['2,3'][-1]['na_mu'], out['pairs']['2,3'][-1]['na_A'], out['pairs']['2,3'][-1]['na_B']), flush=True)
json.dump(out, open(sys.argv[1], 'w'), indent=1)
