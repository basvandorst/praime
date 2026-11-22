"""Round 19 with the exact sums of round 22: where N(30, x) overtakes N(6, x), from N(g, x) = N(1e6; g) +
sum_j (-1)^j U_j(g) * integral_{1e6}^{x} ln(t)^-(j+2) dt (mpmath, 40 digits)."""
import os, json, math
import numpy as np
import mpmath as mp
import common as C
from r22 import U_sums

mp.mp.dps = 40
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
X0 = 10 ** 6
M = C.primes_mask(2, X0 - 2); P = np.append(np.nonzero(M)[0] + 2, 1000003); g0 = np.diff(P)
U = {g: U_sums(g)[0] for g in (6, 30)}


def N(g, lnx):
    s = mp.mpf(int(np.sum(g0 == g)))
    for j in range(len(U[g])):
        s += (-1) ** j * mp.mpf(U[g][j]) * mp.quad(lambda u: mp.e ** u / u ** (j + 2), [mp.log(X0), lnx / 3, lnx])
    return s


u = mp.findroot(lambda v: N(30, v) - N(6, v), mp.mpf(81.1))
x = mp.e ** u
res = {'exact_cross_30_6_log10x': float(u / mp.log(10)), 'x_mantissa_e35': float(x / mp.mpf(10) ** 35)}
print(res)
json.dump(res, open(os.path.join(OUT, 'r19b.json'), 'w'), indent=1)
