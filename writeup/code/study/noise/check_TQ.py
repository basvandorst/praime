"""Finite-H check of the closed forms at the fine modulus Q = 30030 used in (6.1): T_Q(c; H) for c = 0 and c = Q/p * j,
exact sums against (3.1)-(3.2), H = 1e6 (33 multiples of Q), 2e6, 1e8 window sizes are beyond reach; reported as T/H and as
the implied v_Q at lam = ln(6e11)."""
import sys, os, math
sys.argv = ['check_T.py', '1000000']
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import check_T as CT
Q = 30030
Sp = CT.Sprime(Q)
lam = math.log(6e11)
for H in (200000, 1000000):
    t0 = CT.T(Q, 0, H, Sp); c0 = CT.T0_closed(Q, H)
    mu = (Q / 5760) / lam
    print('H=%d c=0: T/H exact %.4f closed %.4f -> v_Q exact %.4f closed %.4f' % (H, t0 / H, c0 / H, 1 - mu + mu * t0 / H, 1 - mu + mu * c0 / H))
    for p in (3, 5, 7, 11, 13):
        c = Q // p
        t = CT.T(Q, c, H, Sp); cl = CT.Tc_closed(Q, c, H)
        print('   c=Q/%d: T/H exact %.5f closed %.5f' % (p, t / H, cl / H))
