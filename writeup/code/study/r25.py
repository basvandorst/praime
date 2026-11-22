"""Round 25 (preregister/round25.txt): round 24 with the polymer overlap correction added to the round-3 law."""
import os, json, math, pickle, time
import numpy as np
from scipy.stats import kstest
import theory2 as T2
import theory_r3 as R3

OUT = R3.OUT
cache = pickle.load(open(R3.CACHE, 'rb'))
P = {int(g): v for g, v in json.load(open(os.path.join(OUT, 'oes_first_occurrence.json')))['P'].items() if (0, int(g)) in cache}
t0 = time.time()
Us = []
for g in sorted(P):
    SB, la = cache[(0, g)]
    pt = T2.Pattern([0, g], 2)                                   # pairs suffice for the correction
    corr = np.array([T2.overlap_correction(pt, l) for l in R3.LAMG])
    u = np.linspace(math.log(1e7), math.log(P[g]), 4001)
    sh = SB / u * np.exp(np.interp(u, R3.LAMG, la + corr))
    Us.append(float(np.trapezoid(sh * np.exp(u) / u, u)))
v = np.array(Us); ks = kstest(v, 'expon')
res = {'g': sorted(P), 'U': Us, 'summary': {'n': len(v), 'mean': float(v.mean()), 'se': float(v.std(ddof=1) / math.sqrt(len(v))),
       'median': float(np.median(v)), 'ks_stat': float(ks.statistic), 'ks_p': float(ks.pvalue)}, 'secs': round(time.time() - t0)}
print(res['summary'], res['secs'])
json.dump(res, open(os.path.join(OUT, 'r25.json'), 'w'), indent=1)
