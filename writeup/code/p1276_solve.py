# Puzzle 1276 (writeup/1276.md): the largest number Q(k) of distinct primes that use every digit 0-9 exactly k times.
# 1. ILP over primes below 10^D, grouped by digit-count vector c: x_c primes are taken from the t_c primes with vector c;
#    sum_c x_c * c = (k, ..., k); maximise sum_c x_c. (HiGHS via scipy, zero gap.)
# 2. Upper bound for all primes by LP duality: for every y >= 0,
#       Q(k) <= k * sum_d y_d + sum_c t_c * max(0, 1 - y.c),
#    and primes with more than D digits add nothing when every such vector has y.c >= 1, which holds if
#       min(y_1, y_3, y_7, y_9) + min(y_1..y_9) + (D - 1) * min(y_0..y_9) >= 1
#    (last digit 1, 3, 7 or 9, first digit nonzero, at least D - 1 further digits). y comes from a dual LP with that
#    condition as a constraint and is rounded to fractions; prob1276 verify recomputes the bound exactly.
# Writes ../data/p1276_solutions.txt: per k a line 'k Q D y0/..,y9' and lines 'c0 .. c9 x' (take any x primes with
# that digit vector, e.g. the smallest).  Run: python3 p1276_solve.py   (about 6 minutes; needs numpy, scipy)
import os, numpy as np
from collections import Counter
from fractions import Fraction as Fr
from scipy.optimize import milp, LinearConstraint, Bounds, linprog

def groups(D):
    N = 10 ** D
    s = np.ones(N, dtype=bool); s[:2] = False
    for i in range(2, int(N ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    cnt = Counter(tuple(str(p).count(d) for d in '0123456789') for p in np.nonzero(s)[0].tolist())
    keys = sorted(cnt)
    return np.array(keys, dtype=float).T, np.array([cnt[c] for c in keys], dtype=float), keys

def dual(C, t, k, D):
    m = C.shape[1]; n = 10 + m + 3                       # variables y (10), z (m), u, v, w
    cost = np.concatenate([np.full(10, k), t, [0, 0, 0]])
    A = [np.hstack([-C.T, -np.eye(m), np.zeros((m, 3))])]; b = [-np.ones(m)]
    for j, ds in enumerate([(1, 3, 7, 9), range(1, 10), range(10)]):
        for d in ds:
            r = np.zeros(n); r[10 + m + j] = 1; r[d] = -1; A.append(r[None]); b.append([0])
    r = np.zeros(n); r[10 + m:] = [-1, -1, -(D - 1)]; A.append(r[None]); b.append([-1])
    res = linprog(cost, A_ub=np.vstack(A), b_ub=np.concatenate(b), bounds=[(0, None)] * n, method='highs')
    return [max(Fr(0), Fr(v).limit_denominator(1000)) for v in res.x[:10]]

KS = {5: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 100, 150, 200],
      6: list(range(250, 1001, 50)),
      7: [2000, 5000, 10000, 20000, 30000, 40000, 50000, 100000]}
out = []
for D, ks in KS.items():
    C, t, keys = groups(D); m = C.shape[1]
    for k in ks:
        res = milp(-np.ones(m), constraints=LinearConstraint(C, np.full(10, k), np.full(10, k)), integrality=np.ones(m),
                   bounds=Bounds(np.zeros(m), t), options={'mip_rel_gap': 0, 'time_limit': 900})
        x = np.round(res.x).astype(int)
        y = dual(C, t, k, D)
        out.append('k %d Q %d D %d y %s' % (k, x.sum(), D, ','.join(str(v) for v in y)))
        out += ['%s %d' % (' '.join(map(str, keys[i])), x[i]) for i in range(m) if x[i]]
        print('k = %d: Q = %d' % (k, x.sum()), flush=True)
here = os.path.dirname(os.path.abspath(__file__))
open(os.path.join(here, '..', 'data', 'p1276_solutions.txt'), 'w').write('\n'.join(out) + '\n')
