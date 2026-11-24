"""Quadrature and truncation checks for k[3](l): GUE at two resolutions; arithmetic part at PB = 1e4 and 1e5;
Q-part numerical against the closed form. usage: kf3_test.py -> kf3_test_out.txt"""
import sys, os, time, math
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
import tplib
LS = [0.5 * i for i in range(1, 17)]
def pick(q, arr): return np.array([arr[int(round(l / q.WS)) - 1] for l in LS])
for kw in (dict(WS=0.25, NS=8, WT=0.5, NT=8),):
    q = tplib.Quad(**kw); g = q.kf3(tplib.T3_gue(0 * q.T, q.T, q.S))
    print('GUE k[3] %s (%d nodes): ' % (kw, len(q.T)) + ' '.join('%.7f' % v for v in pick(q, g)))
q = tplib.Quad()
for L, S in ((9.4, ()), (9.4, (2,))):
    for PB in (10 ** 4,):
        t0 = time.time(); PS = tplib.PrimeSums(PB)
        d = tplib.T3_parts(0 * q.T, q.T, q.S, L, S, PS=PS)
        kQ = pick(q, q.kf3(d['Q'])); kP = pick(q, q.kf3(d['P']))
        print('L=%.1f S=%s PB=%d (%.0fs): Q-part ' % (L, S, PB, time.time() - t0) + ' '.join('%+.6f' % v for v in kQ))
        print('                         P-part ' + ' '.join('%+.6f' % v for v in kP))
    primes = [p for p in tplib.primes_upto(10 ** 6) if p not in S]
    print('                closed-form Q-part ' + ' '.join('%+.6f' % v for v in tplib.kf3_Q_closed(LS, L, primes)))
