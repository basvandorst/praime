"""Round 13 diagnostic 2: spread across residue classes mod 210 of the prime counts n_a and of the pair counts T_a(h)
(p = a, p + h prime and consecutive), each as Var/mean across the admissible classes (1 for Poisson counts)."""
import os, json
import numpy as np
import theory_r3 as R3
import common as C
from math import gcd

Q = 210
out = {}
for win in ('train', 'A', 'B', 'C'):
    L, W = R3.WIN[win]
    M = C.primes_mask(L, W); P = np.nonzero(M)[0].astype(np.int64) + L
    g = np.diff(P); r = P[:-1] % Q
    n_a = np.bincount(r, minlength=Q)[[a for a in range(Q) if gcd(a, Q) == 1]]
    res = {'primes': float(n_a.var(ddof=1) / n_a.mean())}
    for h in (2, 4, 6, 12, 30):
        adm = [a for a in range(Q) if gcd(a, Q) == 1 and gcd(a + h, Q) == 1]
        T = np.bincount(r[g == h], minlength=Q)[adm]
        res['gap%d' % h] = float(T.var(ddof=1) / T.mean())
        res['gap%d_classes' % h] = len(adm)
    out[win] = res
    print(win, {k: round(v, 2) for k, v in res.items()}, flush=True)
json.dump(out, open(os.path.join(R3.OUT, 'r13f.json'), 'w'), indent=1)
