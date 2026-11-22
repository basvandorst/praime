"""Review 26: own count of k-th neighbour repeats in [1e7, 1e8), D <= 200 (plain sieve; spot check with sympy)."""
import math, json
import numpy as np
from sympy import isprime, nextprime
L0, L1 = 10 ** 7, 10 ** 8
s = np.ones(L1, bool); s[:2] = False
for i in range(2, math.isqrt(L1) + 1):
    if s[i]:
        s[i * i::i] = False
P = np.nonzero(s[L0:])[0].astype(np.int64) + L0; del s
assert P[0] == nextprime(L0 - 1) and P[-1] < L1 and nextprime(int(P[-1])) >= L1
rng = np.random.default_rng(1)
for x in rng.integers(0, len(P) - 1, 200):
    assert isprime(int(P[x])) and nextprime(int(P[x])) == int(P[x + 1])
out = {'primes': int(len(P)), 'data': {}}
for k in (1, 2, 3, 4):
    d = P[k:] - P[:-k]; d = d[d <= 200]; out['data'][k] = {'n': int(len(d))}
    for q in (3, 5, 7, 11, 13):
        c = int(np.count_nonzero(d % q == 0)); out['data'][k][q] = c
    # also: hist per D for later
    out['data'][k]['hist'] = np.bincount(d, minlength=201).tolist()
print(len(P), {k: {q: out['data'][k][q] / out['data'][k]['n'] for q in (3, 5, 7, 11, 13)} for k in out['data']})
json.dump(out, open('/home/bas/scratch/study-review26/count.json', 'w'))
