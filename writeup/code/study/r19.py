"""Round 19 (preregister/round19.txt): jumping champions from the cluster law."""
import os, json, math
import numpy as np
from scipy.optimize import brentq
import theory2 as T2

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
P = {6: (T2.Pattern([0, 6], 4), 4), 30: (T2.Pattern([0, 30], 4), 4), 210: (T2.Pattern([0, 210], 3), 3)}


def share(d, u):
    pt, o = P[d]
    return pt.SB / u * math.exp(pt.log_avoid(u, o))


def cum(d, U, U0=math.log(1e3), n=4000):
    """N(d, e^U) up to a constant factor common to all d: integral of e^u share_d(u)/u du (in log scale)."""
    u = np.linspace(U0, U, n)
    v = np.array([share(d, x) for x in u])
    # integrate e^(u - U) to keep numbers finite; the common factor e^U cancels in ratios
    return float(np.trapezoid(np.exp(u - U) * v / u, u))


res = {}
U8 = math.log(1e8)
res['ratio_30_6_at_1e8'] = cum(30, U8) / cum(6, U8)
res['local_ratio_30_6_at_1e8'] = share(30, U8) / share(6, U8)
loc = brentq(lambda u: share(30, u) - share(6, u), 40, 200)
cu = brentq(lambda u: cum(30, u, U0=max(u - 40, math.log(1e3))) - cum(6, u, U0=max(u - 40, math.log(1e3))), 40, 200)
res['local_cross_30_6_log10x'] = loc / math.log(10)
res['cumulative_cross_30_6_log10x'] = cu / math.log(10)
loc2 = brentq(lambda u: share(210, u) - share(30, u), 200, 3000)
cu2 = brentq(lambda u: cum(210, u, U0=u - 60, n=600) - cum(30, u, U0=u - 60, n=600), 200, 3000)
res['local_cross_210_30_log10x'] = loc2 / math.log(10)
res['cumulative_cross_210_30_log10x'] = cu2 / math.log(10)
# the crude model of MERSENNE.md section 52, for comparison: weight S(d) e^{-d/u}
S = {d: T2.S1([0, d]) for d in P}
crude = lambda d, u: S[d] * math.exp(-d / u)
res['crude_local_cross_30_6_log10x'] = brentq(lambda u: crude(30, u) - crude(6, u), 10, 500) / math.log(10)
res['crude_local_cross_210_30_log10x'] = brentq(lambda u: crude(210, u) - crude(30, u), 100, 5000) / math.log(10)
for k, v in res.items():
    print(k, round(v, 4))
json.dump(res, open(os.path.join(OUT, 'r19.json'), 'w'), indent=1)
