"""Theory only: scatter of the swapped-term predictions from the Euler-product cutoff of H (5e4, as frozen, against 2e5
and 1e4) and from the u-step (0.06 against 0.03), on 12 near pairs of the fresh design NN. -> swapcheck_out.txt"""
import os, sys, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from crosslib import *
from swaplib import cswap
import swap_run
W = 0.3 * np.arange(1, 17)
PL = swap_run.pairs('NN'); PL = sorted(PL, key=lambda x: -x[2])
sel = PL[::max(1, len(PL) // 12)][:12]
res = {}
for PM in (10000, 50000, 200000):
    P = primes_upto(PM)
    for du in ((0.06, 0.03) if PM == 50000 else (0.06,)):
        out = []
        for a, b, f in sel:
            x = kron_table([a, b], P); out.append(cswap(a, b, W, 100.0, 600.0, 0.02, P, x[0], x[1], du=du))
        res[(PM, du)] = np.array(out)
ref = res[(200000, 0.06)]
for key, v in res.items():
    d = v - ref; rel = np.sqrt((d ** 2).sum() / (ref ** 2).sum())
    print('P<=%6d du=%.2f: rms difference from (2e5, 0.06) %.5f, relative %.3f; slope on it %.3f' % (key[0], key[1], np.sqrt((d ** 2).mean()), rel, float((v * ref).sum() / (ref * ref).sum())))
print('pairs: ' + ', '.join('(%d,%d,F=%.2f)' % s for s in sel))
print('reference values at w = 0.3, 1.5, 3.0, 4.5 (rows = pairs):')
for s, r in zip(sel, ref): print('  (%d,%d) ' % (s[0], s[1]) + ' '.join('%+.5f' % v for v in r[[0, 4, 9, 14]]))
