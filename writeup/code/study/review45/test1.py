import numpy as np, sys
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
from mylaw import *
Sf = lambda r: np.sinc(r)
rng = np.random.default_rng(1)
t = rng.uniform(0.01, 4, 20); s = t + rng.uniform(0.01, 4, 20)
for L in (9.5, 11.0):
    g = T3_P(t, s, L, [], 1.0, gue=True)
    print('GUE skeleton max diff', np.max(np.abs(g - 2*Sf(t)*Sf(s)*Sf(s - t))))
# coincidence check: T3(0,t,s) -> -2 T2(s) as t->0  (full = P + Q)
for D in (1009, 1004, -1155, 1001):
    bad = badset(D); kp2 = kphi2(D); L = 10.5
    ss = np.array([0.7, 1.9, 3.3, 5.6])
    for tt in (1e-3, 1e-4):
        tv = np.full(4, tt)
        f = T3_P(tv, ss, L, bad, kp2, P=20000) + T3_Q(tv, ss, L, bad)
        print(D, tt, np.round(f + 2*T2(ss, L, bad, kp2), 6))
