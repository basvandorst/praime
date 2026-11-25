"""check_vlib.py: vlib (extended tables) against the repo's bklib/tplib (registered law) on their common range eps <= 7.6,
and the exact constraint T3(0, 0, v) = -2 T2(v) (R3 = 0 at a coincidence) far beyond it."""
import sys, os, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "lib")); sys.path.insert(0, HERE)
import tplib, bklib, vlib
rng = np.random.default_rng(1)
L = 10.0
r = rng.uniform(0.01, 7.5 * L / (2 * math.pi), 200)
for S in [(), (2,), (3, 5)]:
    a = tplib.bklib.C_of_r(tplib.G, L, S, r); b = vlib.C_of_r(L, S, r)
    print('C_of_r  S=%-7s max |vlib - bklib| = %.2e  (max |C| %.2e)' % (S, np.abs(a - b).max(), np.abs(a).max()))
PS_t = tplib.PrimeSums(2000); PS_v = vlib.PrimeSums(2000)
n = 300; u = np.zeros(n); v = rng.uniform(-5, 5, n); w = rng.uniform(-5, 5, n)
for S in [(), (2,), (3,), (2, 7)]:
    a = tplib.T3_parts(u, v, w, L, S, PS_t)['P']
    B1 = [PS_v.B1_all(x, y) for (x, y) in vlib.b1_args(u, v, w, L)]
    b = vlib.T3P(u, v, w, L, [S], B1)[0]
    print('T3 P-part S=%-7s max |vlib - tplib| = %.2e  (max |P| %.2e)' % (S, np.abs(a - b).max(), np.abs(a).max()))
# coincidence constraint far out: T3(0, t, v) -> -2 T2(v) as t -> 0, with T2 = T2_GUE + C
for vv in (20.0, 45.0, 90.0):
    for S in [(), (2,)]:
        t = 1e-4; uu = np.array([0.0]); x = np.array([t]); y = np.array([vv])
        B1 = [PS_v.B1_all(a1, b1) for (a1, b1) in vlib.b1_args(uu, x, y, L)]
        P = vlib.T3P(uu, x, y, L, [S], B1)[0][0]
        Qp = 0.0
        k = 2 * math.pi / L; pos = {'u': 0.0, 'v': t * k, 'w': vv * k}
        Qall = 0
        for (A, B1n, B2n) in vlib.TERMS:
            al = -pos[A]; e1 = al + pos[B1n]; e2 = al + pos[B2n]
            Qall += PS_v_Q if False else sum(-(lp ** 3) / ((p * np.exp(1j * e1 * lp) - 1) * (p * np.exp(1j * e2 * lp) - 1))
                                             for p, lp in zip(PS_v.P, PS_v.lnP) if p not in S)
        Q = (2 / L ** 3 * Qall).real
        G = vlib.T3_gue(0.0, t, vv)
        T2 = vlib.T2_gue(vv) + vlib.C_of_r(L, S, np.array([vv]))[0]
        print('v=%5.1f S=%-5s T3(0,1e-4,v) = %+.4e   -2 T2(v) = %+.4e' % (vv, S, G + Q + P, -2 * T2))
