"""Round 23 (preregister/round23.txt): the exact Hardy-Littlewood gap law against Oliveira e Silva's counts of each
gap below 4e18. N_law = N_true(1e6; g) + integral_{1e6}^{4e18} sum_j (-1)^j U_j(g) / ln(t)^(j+2) dt."""
import os, json, math
import numpy as np
import mpmath as mp
import common as C
from r22 import U_sums

mp.mp.dps = 30
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
X0, X = 10 ** 6, 4 * 10 ** 18
data = json.load(open(os.path.join(OUT, 'oes_gap_counts_4e18.json')))['N']

M = C.primes_mask(2, X0 - 2); P = np.nonzero(M)[0] + 2
P = np.append(P, 1000003)                                   # the prime after 1e6 closes the last gap below 1e6
g0 = np.diff(P)
Ilog = {}


def I(k):
    if k not in Ilog:
        Ilog[k] = mp.quad(lambda u: mp.e ** u / u ** k, [mp.log(X0), mp.log(1e9), mp.log(1e14), mp.log(X)])
    return Ilog[k]


out = {}
for g in range(2, 72, 2):
    r = U_sums(g)
    if not isinstance(r, tuple):
        continue
    U, _ = r
    law = sum(((-1) ** j) * mp.mpf(U[j]) * I(j + 2) for j in range(len(U))) + int(np.sum(g0 == g))
    d = data[str(g)]
    diff = d - float(law)
    out[g] = {'data': d, 'law': float(law), 'rel': diff / d, 'z_sqrtN': diff / math.sqrt(d),
              'z_oes': diff / math.sqrt(2 * d * math.log(math.log(d)))}
    print(g, d, '%.0f' % float(law), '%+.2e' % (diff / d), 'z %.1f' % out[g]['z_sqrtN'], 'oes %.2f' % out[g]['z_oes'], flush=True)
json.dump(out, open(os.path.join(OUT, 'r23.json'), 'w'), indent=1)
