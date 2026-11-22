"""Check 6: consecutive twin primes. Exact HL inclusion-exclusion over all subsets of the live interior pair events
for B = {0, 2, h, h+2}, h <= 120, versus the table values (order 1 / order 3) at train.
share(lam) = P(B) / (S(0,2)/lam^2), P(B) = sum_A (-1)^|A| S(B u (uA)) / lam^(4 + |uA \\ B|);
window average with weight 1/ln^2 t over [1e7, 1e8)."""
import math, itertools
import numpy as np
from myS import S

L, W = 10**7, 9 * 10**7
TABLE = {6: (0.010146, 0.010146), 12: (0.027055, 0.027055), 18: (0.019746, 0.019743), 24: (0.012319, 0.012319),
         30: (0.038235, 0.038217), 36: (0.010867, 0.010865), 42: (0.034844, 0.034807), 48: (0.018627, 0.018585),
         60: (0.027358, 0.027285), 90: (0.027141, 0.026927), 120: (0.021143, 0.020902)}
MEAS = {6: 0.010146, 12: 0.027097, 18: 0.019833, 24: 0.012265, 30: 0.038088, 36: 0.010998, 42: 0.034584,
        48: 0.018695, 60: 0.027417, 90: 0.026958, 120: 0.020667}


def coeffs(h):
    B = [0, 2, h, h + 2]
    SB = S(B)
    ev = [m for m in range(1, h) if S(sorted(set(B) | {m, m + 2}), P=h + 2) > 0]
    # c[j] = coefficient of lam^-j (j = number of new elements)
    c = {}
    n = len(ev)
    # enumerate subsets with pruning: once dead, supersets dead
    def rec(i, cur, sign):
        val = S(sorted(cur), P=h + 2)
        if val == 0:
            return
        j = len(cur) - 4
        c[j] = c.get(j, 0.0) + sign * val / SB
        for k in range(i, n):
            m = ev[k]
            rec(k + 1, cur | {m, m + 2}, -sign)
    rec(0, set(B), 1)
    return SB, ev, c


def share(h):
    SB, ev, c = coeffs(h)
    t = np.linspace(L, L + W, 200001)
    lam = np.log(t)
    E = sum(v * lam ** (-j) for j, v in c.items())
    f = SB / S([0, 2]) / lam ** 2 * E
    wgt = 1 / lam ** 2
    # order 1 (independent interior events): prod (1 - rho_e / lam^2)
    o1 = np.ones_like(lam)
    for m in ev:
        o1 *= 1 - S(sorted({0, 2, h, h + 2, m, m + 2}), P=h + 2) / SB / lam ** 2
    f1 = SB / S([0, 2]) / lam ** 2 * o1
    avg = lambda g: float(np.trapezoid(g * wgt, t) / np.trapezoid(wgt, t))
    return len(ev), avg(f), avg(f1)


if __name__ == '__main__':
    print(' h  live  exact HL     my order1   table o1 / o3          measured   exact/o3-1')
    for h in [x for x in TABLE if x <= 90]:
        n, ex, o1 = share(h)
        print('%3d  %3d  %.6f   %.6f   %.6f / %.6f   %.6f   %+.2e' % (h, n, ex, o1, *TABLE[h], MEAS[h], ex / TABLE[h][1] - 1),
              flush=True)
