"""The prediction: kappa3(l) = k[3](l) + 3 Sigma2(l) - 2 l of the window count of Dirichlet L-zeros (ratios conjecture /
Bogomolny-Keating, all arithmetic terms), its GUE value, and the per-prime change Delta kappa3(l; p, L, S) with its parts.
usage: kappa3_tables.py -> kappa3_tables_out.txt, kappa3_tables.npz"""
import sys, os, math, json
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
import tplib, bklib
q = tplib.Quad(); PS = tplib.PrimeSums(10 ** 4)
LG = q.lgrid; LS = [0.5 * i for i in range(1, 17)]; IX = [int(round(l / q.WS)) - 1 for l in LS]
PRIMES = tplib.primes_upto(10 ** 6)
kg = q.kf3(tplib.T3_gue(0 * q.T, q.T, q.S)); s2g = tplib.sigma2_gue(LG); k3g = kg + 3 * s2g - 2 * LG
f = lambda a: ' '.join('%+.5f' % v for v in np.asarray(a)[IX])
print('l           ' + ' '.join('%8.1f' % l for l in LS))
print('GUE k[3]    ' + f(kg)); print('GUE Sigma2  ' + f(s2g)); print('GUE kappa3  ' + f(k3g))
# derivative of Sigma2_GUE on the fine grid (for the coupling heuristic)
s2g_fine = bklib.sigma2_gue(); ds2g = np.gradient(s2g_fine, bklib.LGRID); ds2g = ds2g[[int(round(l / bklib.W)) - 1 for l in LG]]
out = {'l': LG.tolist(), 'gue_kf3': kg.tolist(), 'gue_s2': s2g.tolist(), 'gue_k3': k3g.tolist(), 'cases': {}}
def full(L, S):
    d = tplib.T3_parts(0 * q.T, q.T, q.S, L, S, PS=PS)
    kQ = tplib.kf3_Q_closed(LG, L, [p for p in PRIMES if p not in S]); kP = q.kf3(d['P'])
    kf = kg + kQ + kP; s2 = tplib.sigma2_bk(L, S, LG)
    return dict(kf3=kf, kQ=kQ, kP=kP, s2=s2, k3=kf + 3 * s2 - 2 * LG)
for L in (8.1, 9.4, 10.2, 10.7):
    base = full(L, ())
    print('\n=== L = %.1f ===' % L)
    print('k = 1 (no bad prime): kappa3        ' + f(base['k3']))
    print('                      kappa3 - GUE  ' + f(base['k3'] - k3g))
    print('   parts: 3(Sigma2 - GUE)           ' + f(3 * (base['s2'] - s2g)))
    print('          k[3] Q-part (diagonal)    ' + f(base['kQ']))
    print('          k[3] P-part minus GUE     ' + f(base['kP']))
    out['cases']['%.1f|' % L] = {k: v.tolist() for k, v in base.items()}
    for p in (2, 3, 5, 7):
        c = full(L, (p,)); out['cases']['%.1f|%d' % (L, p)] = {k: v.tolist() for k, v in c.items()}
        dk3 = c['k3'] - base['k3']; ds2 = c['s2'] - base['s2']
        ds2d = np.array([bklib.dS_diag(p, L, l) for l in LG]); ds2o = ds2 - ds2d
        dQ = c['kQ'] - base['kQ']; dP = c['kP'] - base['kP']
        print('p = %d: Delta kappa3                 ' % p + f(dk3))
        print('       3 Delta Sigma2 (diag+off)     ' + f(3 * ds2) + '   [diag ' + f(3 * ds2d) + ']')
        print('       Delta k[3] Q (closed form)    ' + f(dQ))
        print('       Delta k[3] P                  ' + f(dP))
        print('       3 Delta Sigma2 + Delta k[3]P  ' + f(3 * ds2 + dP))
        print('       coupling 3 Sigma2_GUE\'(l) dS2d ' + f(3 * ds2g * ds2d))
json.dump(out, open('kappa3_tables.json', 'w'))
