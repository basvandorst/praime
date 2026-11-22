"""Round 24 (preregister/round24.txt): first occurrences of large gaps against the cluster law.
U_g = N_law(P(g); g) = integral from 1e7 to P(g) of share_g(u) e^u / u du (no gap g >= 180 occurs below 1e7)."""
import os, json, math, pickle
import numpy as np
from scipy.stats import kstest
import theory2 as T2
import theory_r3 as R3

OUT = R3.OUT
cache = pickle.load(open(R3.CACHE, 'rb'))
P = {int(g): v for g, v in json.load(open(os.path.join(OUT, 'oes_first_occurrence.json')))['P'].items() if (0, int(g)) in cache}   # the cache ends at g = 454


def U(g, x, crude=False):
    u = np.linspace(math.log(1e7), math.log(x), 4001)
    if crude:
        sh = T2.S1([0, g]) * np.exp(-g / u) / u
    else:
        SB, la = cache[(0, g)]
        sh = SB / u * np.exp(np.interp(u, R3.LAMG, la))
    return float(np.trapezoid(sh * np.exp(u) / u, u))


res = {'g': [], 'P': [], 'U': [], 'U_crude': []}
for g in sorted(P):
    res['g'].append(g); res['P'].append(P[g]); res['U'].append(U(g, P[g])); res['U_crude'].append(U(g, P[g], True))
for k in ('U', 'U_crude'):
    v = np.array(res[k])
    ks = kstest(v, 'expon')
    res[k + '_summary'] = {'n': len(v), 'mean': float(v.mean()), 'se': float(v.std(ddof=1) / math.sqrt(len(v))),
                           'median': float(np.median(v)), 'ks_stat': float(ks.statistic), 'ks_p': float(ks.pvalue)}
    print(k, res[k + '_summary'])
json.dump(res, open(os.path.join(OUT, 'r24.json'), 'w'), indent=1)
