"""Independent reference for prime9.html (plain Python integers, sympy for primality): for the primes up to X and the
Mersenne exponents, the depth of Gilbreath's waterfall (first row d >= 1 whose first entry is 1 and all later entries
are 0 or 2), whether the first column is 1 down to there, and the first column. Also checks the triangle against
OEIS A036277 (position of the first entry > 2 in row n, rows 0..50).   usage: ref9.py out.json"""
import sys, json
from sympy import isprime, primerange
A036277 = [2,4,9,15,15,26,25,24,23,26,60,99,98,99,98,175,177,177,177,177,292,291,290,741,875,874,873,874,873,872,871,870,869,868,867,2181,2180,2179,2178,2772,2771,2770,2769,2768,2767,2766,2765,2764,2764,2764,2764,3367]
def run(seq):
    row = list(seq); d = 0; firsts = [row[0]]
    while True:
        if d > 0:
            if row[0] != 1: return dict(depth=d, holds=False, firsts=firsts)
            if len(row) <= 1 or all(v in (0, 2) for v in row[1:]): return dict(depth=d, holds=True, firsts=firsts)
        if len(row) <= 1: return dict(depth=d, holds=True, firsts=firsts)
        row = [abs(row[i + 1] - row[i]) for i in range(len(row) - 1)]; d += 1; firsts.append(row[0])
out = {}
for X in (1000, 5000, 100000, 1000000):
    P = list(primerange(2, X + 1)); r = run(P); out['primes_%d' % X] = dict(depth=r['depth'], holds=r['holds'], firsts=r['firsts'][:12], n=len(P))
M = [2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243,110503,132049,216091,756839,859433]
assert all(isprime(p) for p in M)
r = run(M); out['mexp_1000000'] = dict(depth=r['depth'], holds=r['holds'], firsts=r['firsts'][:12], n=len(M))
row = list(primerange(2, 10 ** 5)); pos = []
for n in range(len(A036277)):
    pos.append(next(i for i, v in enumerate(row) if v > 2) + 1)
    row = [abs(row[i + 1] - row[i]) for i in range(len(row) - 1)]
out['A036277_rows_0_to_51_agree'] = pos == A036277
json.dump(out, open(sys.argv[1], 'w')); print(json.dumps(out))
