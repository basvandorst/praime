# Puzzle 1276 (writeup/1276.md): independent check of ../data/p1276_solutions.txt.
# For every k: rebuild the set of primes (the x smallest primes with each listed digit vector), check that they are
# distinct primes using every digit 0-9 exactly k times, and that there are Q of them. Then recompute from scratch,
# in exact fractions, the LP-duality bound valid for all primes:
#     Q(k) <= k * sum(y) + sum over digit vectors c of t_c * max(0, 1 - y.c)    (t_c = number of primes < 10^D with vector c)
# and check that primes with more than D digits cannot add to it:
#     min(y_1, y_3, y_7, y_9) + min(y_1..y_9) + (D - 1) min(y) >= 1.
# Q = floor(bound) proves Q optimal. Run: python3 p1276_verify.py   (about 1 minute; needs sympy)
import os, random
from fractions import Fraction as Fr
from collections import defaultdict
from sympy import isprime

def primes_below(N):
    s = bytearray([1]) * N; s[0:2] = b'\x00\x00'
    for i in range(2, int(N ** 0.5) + 1):
        if s[i]: s[i * i::i] = bytearray(len(range(i * i, N, i)))
    return [i for i in range(N) if s[i]]

here = os.path.dirname(os.path.abspath(__file__))
blocks = []
for line in open(os.path.join(here, '..', 'data', 'p1276_solutions.txt')):
    w = line.split()
    if w[0] == 'k':
        blocks.append({'k': int(w[1]), 'Q': int(w[3]), 'D': int(w[5]), 'y': [Fr(v) for v in w[7].split(',')], 'x': []})
    else:
        blocks[-1]['x'].append((tuple(map(int, w[:10])), int(w[10])))

page = {1: 6, 2: 10, 3: 14, 4: 18, 5: 22, 6: 26, 7: 29, 8: 33, 9: 36, 10: 39, 20: 73, 30: 103, 40: 131, 50: 159,
        100: 295, 150: 423, 200: 543, 250: 659, 300: 776, 350: 893, 400: 1008, 450: 1117, 500: 1226, 550: 1335,
        600: 1444, 650: 1553, 700: 1662, 750: 1769, 800: 1875, 850: 1980, 900: 2084, 950: 2185, 1000: 2285,
        2000: 4285, 5000: 9895, 10000: 18503, 20000: 35168, 30000: 51515, 40000: 67251, 50000: 82709, 100000: 155645}
cache = {}
for b in blocks:
    k, Q, D, y = b['k'], b['Q'], b['D'], b['y']
    if D not in cache:
        byvec = defaultdict(list)
        for p in primes_below(10 ** D):
            byvec[tuple(str(p).count(d) for d in '0123456789')].append(p)
        cache[D] = byvec
    byvec = cache[D]
    S = [p for c, x in b['x'] for p in byvec[c][:x]]
    assert all(len(byvec[c]) >= x for c, x in b['x'])
    assert len(S) == len(set(S)) == Q
    assert all(isprime(p) for p in random.sample(S, min(len(S), 2000)))     # spot check of the sieve
    tot = [0] * 10
    for p in S:
        for ch in str(p): tot[int(ch)] += 1
    assert tot == [k] * 10, (k, tot)
    assert all(v >= 0 for v in y)
    assert min(y[d] for d in (1, 3, 7, 9)) + min(y[1:]) + (D - 1) * min(y) >= 1
    B = k * sum(y)
    for c, ps in byvec.items():
        s = sum(y[d] * c[d] for d in range(10))
        if s < 1: B += len(ps) * (1 - s)
    assert B // 1 == Q, (k, Q, float(B))
    tag = 'equals the page' if page[k] == Q else 'page has %d' % page[k]
    print('k = %6d: %6d primes, every digit used %d times; bound %.4f, so %d is optimal (%s)' % (k, Q, k, float(B), Q, tag))
print('all checks passed')
