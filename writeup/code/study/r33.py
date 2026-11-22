"""Round 33 (diagnostic, no data): the order-5 term against the order-4 term for the long gap 180 at lam = 18, 20, 24,
to see how fast the Bahadur series falls for long gaps (r8 had order 5 only for spans <= 120)."""
import os, json, math, time
from itertools import combinations, permutations
import numpy as np
import r8

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
LAMS = [18.0, 20.0, 24.0]
g = 180
t0 = time.time()
SB, live, s1, s2, s3, Ssub, n = r8.tables([0, g])
Q4 = np.array(list(combinations(range(n), 4)), dtype=np.int64)
s4v = np.concatenate([Ssub(live[Q4[i:i + 500000]]) for i in range(0, len(Q4), 500000)])
t4 = np.zeros(len(LAMS))
a, b, c, d = Q4.T
K4 = (s4v - s3[a, b, c] * s1[d] - s3[a, b, d] * s1[c] - s3[a, c, d] * s1[b] - s3[b, c, d] * s1[a]
      - s2[a, b] * s2[c, d] - s2[a, c] * s2[b, d] - s2[a, d] * s2[b, c]
      + 2 * (s2[a, b] * s1[c] * s1[d] + s2[a, c] * s1[b] * s1[d] + s2[a, d] * s1[b] * s1[c]
             + s2[b, c] * s1[a] * s1[d] + s2[b, d] * s1[a] * s1[c] + s2[c, d] * s1[a] * s1[b])
      - 6 * s1[a] * s1[b] * s1[c] * s1[d])
for li, lam in enumerate(LAMS):
    q = 1 - s1 / lam
    t4[li] = float(np.sum(K4 / lam ** 4 / (q[a] * q[b] * q[c] * q[d])))
del K4
t5 = -r8.order5_term(live, s1, s2, s3, s4v, Q4, Ssub, n, np.array(LAMS))      # r8 returns the signed term; keep its sign
res = {'g': g, 'live': int(n), 'order4_term_ln': dict(zip(map(str, LAMS), t4.tolist())),
       'order5_term_ln': dict(zip(map(str, LAMS), (-t5).tolist())), 'secs': round(time.time() - t0)}
print(res, flush=True)
json.dump(res, open(os.path.join(OUT, 'r33.json'), 'w'), indent=1)
