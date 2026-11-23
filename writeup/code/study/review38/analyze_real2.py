import numpy as np
from common import *
z = np.load('/home/bas/scratch/review-logic3/real.npz')
Ds, ls, L = z['Ds'], z['ls'], z['L']
n = len(Ds)
I = np.array([[1.0 if abs(D) % p == 0 else 0.0 for p in PR] for D in Ds])
X = np.column_stack([np.ones(n), I])
# unmodeled primes > 13
for key in ['dense', 'rand', 'trim']:
    S = z[key]
    for j, l in enumerate(ls):
        if l not in (1.0, 2.0, 3.5, 6.5, 7.0, 7.5): continue
        rest = np.array([sum(dS(l, p, L[k]) for p in primes_of(D) if p > 13) for k, D in enumerate(Ds)])
        b0, se0 = ols(X, S[:, j])
        b1, se1 = ols(np.column_stack([X, rest]), S[:, j])
        b2, _ = ols(X, S[:, j] - rest)  # offset with coefficient 1
        print(key, l, 'base', np.round(b0[1:], 4), '| +rest', np.round(b1[1:7], 4), 'coef rest %.2f+-%.2f' % (b1[7], se1[7]),
              '| offset', np.round(b2[1:], 4))
# corr between indicators and rest at l=1
rest = np.array([sum(dS(1, p, L[k]) for p in primes_of(D) if p > 13) for k, D in enumerate(Ds)])
print('mean rest by 2|D', rest[I[:, 0] == 1].mean(), rest[I[:, 0] == 0].mean())
for i, p in enumerate(PR):
    print(p, 'n=%d' % I[:, i].sum(), 'mean rest in group %.4f out %.4f' % (rest[I[:, i] == 1].mean(), rest[I[:, i] == 0].mean()))
# half-integer vs integer alternation, check full 96 slope split
