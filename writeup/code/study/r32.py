"""Round 32 (diagnostic, no data): size of the order-4 Bahadur term for long gaps at small lam, chunked over the
quadruples, for gaps 240, 242, 300, 302 at lam = 18, 20, 24 (the sizes of the first occurrences of rounds 24-25)."""
import os, sys, json, math, time
from itertools import combinations
import numpy as np
import r8

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
LAMS = [18.0, 20.0, 24.0]
path = os.path.join(OUT, 'r32.json')
res = json.load(open(path)) if os.path.exists(path) else {}
for g in [int(x) for x in sys.argv[1:]]:
    if str(g) in res:
        continue
    t0 = time.time()
    SB, live, s1, s2, s3, Ssub, n = r8.tables([0, g])
    tot = np.zeros(len(LAMS))
    it = combinations(range(n), 4)
    while True:
        Q = np.array([x for _, x in zip(range(500000), it)], dtype=np.int64)
        if len(Q) == 0:
            break
        s4v = Ssub(live[Q]); a, b, c, d = Q.T
        K4 = (s4v - s3[a, b, c] * s1[d] - s3[a, b, d] * s1[c] - s3[a, c, d] * s1[b] - s3[b, c, d] * s1[a]
              - s2[a, b] * s2[c, d] - s2[a, c] * s2[b, d] - s2[a, d] * s2[b, c]
              + 2 * (s2[a, b] * s1[c] * s1[d] + s2[a, c] * s1[b] * s1[d] + s2[a, d] * s1[b] * s1[c]
                     + s2[b, c] * s1[a] * s1[d] + s2[b, d] * s1[a] * s1[c] + s2[c, d] * s1[a] * s1[b])
              - 6 * s1[a] * s1[b] * s1[c] * s1[d])
        for li, lam in enumerate(LAMS):
            q = 1 - s1 / lam
            tot[li] += float(np.sum(K4 / lam ** 4 / (q[a] * q[b] * q[c] * q[d])))
    res[str(g)] = {'live': int(n), 'order4_term_ln': dict(zip([str(l) for l in LAMS], tot.tolist())), 'secs': round(time.time() - t0)}
    json.dump(res, open(path, 'w'), indent=1)
    print(g, res[str(g)], flush=True)
