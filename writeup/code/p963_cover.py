# Puzzle 963 (writeup/963.md): fewest arithmetic progressions of primes (>= 3 terms, all terms among the first X odd
# primes) covering the first X odd primes. Set cover over the MAXIMAL progressions (HiGHS via scipy).
# Usage: python3 p963_cover.py X [time_limit_s] [loose]   (loose: progressions may continue past the X-th odd prime)
import sys, numpy as np
from sympy import prime
from scipy.optimize import milp, LinearConstraint, Bounds, linprog
from scipy.sparse import lil_matrix

X = int(sys.argv[1]); tl = float(sys.argv[2]) if len(sys.argv) > 2 else 600
P = [prime(i) for i in range(2, X + 2)]          # first X odd primes
S = set(P); idx = {p: i for i, p in enumerate(P)}
aps = []
for i, a in enumerate(P):
    for b in P[i + 1:]:
        d = b - a
        if a - d in S: continue                    # not maximal to the left
        if b + d not in S: continue                # fewer than 3 terms
        t = [a, b]
        while t[-1] + d in S: t.append(t[-1] + d)
        aps.append(t)
from sympy import isprime
if 'loose' in sys.argv:          # progressions may continue past the first X odd primes: adds pairs {a, b} with 2b - a prime
    for i, a in enumerate(P):
        for b in P[i + 1:]:
            if a - (b - a) not in S and b + (b - a) not in S and isprime(2 * b - a): aps.append([a, b])
print('X =', X, 'max prime', P[-1], 'covering sets:', len(aps))
m = len(aps)
A = lil_matrix((X, m))
for j, t in enumerate(aps):
    for p in t: A[idx[p], j] = 1
A = A.tocsr()
lp = linprog(np.ones(m), A_ub=-A, b_ub=-np.ones(X), bounds=(0, 1), method='highs')
print('LP bound %.3f' % lp.fun)
res = milp(np.ones(m), constraints=LinearConstraint(A, 1, np.inf), integrality=np.ones(m), bounds=Bounds(0, 1),
           options={'time_limit': tl, 'mip_rel_gap': 0})
print(res.message, ' cover size', round(res.fun), ' lower bound', getattr(res, 'mip_dual_bound', None))
sol = [aps[j] for j in range(m) if res.x[j] > 0.5]
assert set(p for t in sol for p in t) == S
open('p963_cover_%d.out' % X, 'w').write('\n'.join(' '.join(map(str, t)) for t in sol) + '\n')
