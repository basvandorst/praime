"""Item 4: is v_2310 = v_30030 * (1 + 10 * mean eps over pairs of 30030-subclasses inside one 2310-class)?
(class A mod 2310 holds 11 admissible subclasses mod 30030; X_A = sum of them, mean 11 mu.)"""
import sys, math, json
import numpy as np
import theory as T

L0 = int(float(sys.argv[1])); W = int(float(sys.argv[2]))
Q = 30030
S = T.setup(Q, L0, W)
r = S['r']; mu = float(r.sum())
F0, _, _ = T.F(S, 0)
var = float(np.sum(r * (1 - r))) + F0
v30030 = var / mu
A = 17                                              # admissible class mod 2310 (17, 19 units)
subs = [a for a in range(A, Q, 2310) if math.gcd(a * (a + 2), Q) == 1]
assert len(subs) == 11
Fc = {}
cov = 0.0
for a in subs:
    for b in subs:
        if a == b: continue
        c = (b - a) % Q
        if c not in Fc:
            Fc[c] = T.F(S, c)[0]
        cov += Fc[c]
v2310_from = (11 * var + cov) / (11 * mu)
S2 = T.setup(2310, L0, W)
r2 = S2['r']; F2, _, _ = T.F(S2, 0)
v2310 = (float(np.sum(r2 * (1 - r2))) + F2) / float(r2.sum())
# eps for pairs in the same class mod 30030/p, other p
other = {}
for p in (5, 7, 11, 13):
    c = Q // p
    other[p] = T.F(S, c)[0] / var if c % 6 == 0 else None
out = dict(L0=L0, W=W, v30030=v30030, v2310_direct=v2310, v2310_from_30030=v2310_from,
           rel_diff=v2310_from / v2310 - 1, mean_eps_within=cov / (110 * var), eps_by_p=other)
print(json.dumps(out, indent=1))
json.dump(out, open(f'consist_{L0}.json', 'w'), indent=1)
