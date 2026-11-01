# Puzzle 207 Q2: every inventory sequence ends in a cycle of period 1, 2, 3 or 4; the largest period is 4
# (writeup/207.md). Independent of the C programs p207_phi.c / p207_orders.c.
# Run: python3 p207_verify.py          (count cycles, order dynamics for supports <= 8, the data files; ~30 s)
#      python3 p207_verify.py --full   (all supports, and the count of all self-inventoried numbers; ~20 min)
import sys
from functools import lru_cache
from itertools import permutations, combinations, combinations_with_replacement

def g(s):
    """Inventory step on a digit string: for each digit in order of first appearance, its count then the digit."""
    order = []
    for ch in s:
        if ch not in order: order.append(ch)
    return ''.join(str(s.count(d)) + d for d in order)

# 0. The page's examples.
s = '1'
for _ in range(15): s = g(s)
assert s in ('23322114', '32232114') and g(g(s)) == s
assert g('22') == '22'
x = '40'
seq = [x]
for _ in range(20): x = g(x); seq.append(x)
assert seq[-1] == seq[-5] and len(set(seq[-4:])) == 4                       # 40 ends in a cycle of period 4

# 1. Length bound: the image of a number with L digits has at most F(L) digits (exact maximum over count vectors).
@lru_cache(None)
def best(L, k):                       # max sum of (1 + len(count)) over at most k present digits with counts summing to L
    if L == 0: return 0
    if k == 0: return -10**9
    return max(best(L, k - 1), max(1 + len(str(c)) + best(L - c, k - 1) for c in range(1, L + 1)))
F = {L: best(L, 10) for L in range(1, 60)}
assert all(F[L] <= 21 for L in range(1, 22)) and all(F[L] < L for L in range(22, 60))
assert F[19] == 21                                       # e.g. counts (10, 1, ..., 1): 19 digits -> 21 digits
# For L >= 60, F(L) <= 10 + 10 * len(str(L)) < L, so every sequence eventually stays at <= 21 digits.
assert all(10 + 10 * len(str(L)) < L for L in range(60, 10**4))

# 2. Count dynamics: the digit counts of g(n) depend only on the counts c of n:
#    Phi(c)[e] = [c_e > 0] + (number of digits e in the decimal counts c_d, d present).
def phi(c):
    o = [0] * 10
    for d in range(10):
        if c[d]:
            o[d] += 1
            for ch in str(c[d]): o[int(ch)] += 1
    return tuple(o)
# Every vector in a cycle is an image of Phi; an image is determined by the support S and the multiset of counts.
images = set()
for k in range(1, 11):
    parts = [m for m in combinations_with_replacement(range(1, 22), k) if sum(m) <= 21]
    for S in combinations(range(10), k):
        for m in parts:
            v = [0] * 10
            for d in S: v[d] += 1
            for cnt in m:
                for ch in str(cnt): v[int(ch)] += 1
            if sum(v) <= 21: images.add(tuple(v))
images.add(tuple([0] * 10))
cycles, seen = [], {}
for v0 in images:
    path, x = [], v0
    while x not in seen and x not in path and sum(x) <= 21:
        path.append(x); x = phi(x)
    if x in path:
        cyc = path[path.index(x):]
        key = min(cyc)
        if key not in seen: cycles.append(cyc)
    for y in path: seen[y] = True
cycles = {min(c): c for c in cycles}
assert len(cycles) == 151 and max(len(c) for c in cycles.values()) == 3
print('count map: %d images with digit sum <= 21, %d cycles, count periods %s'
      % (len(images), len(cycles), sorted({len(c) for c in cycles.values()})))

# 3. Order dynamics: the state (order of first appearance, counts) determines the next number. For each count cycle
#    and every order of the support of its first vector, follow one turn and find the cycle lengths of that map.
def turn(order, cyc):
    for c in cyc:
        s = ''.join(str(c[d]) + str(d) for d in order)
        nc = [0] * 10
        o = []
        for ch in s:
            nc[int(ch)] += 1
            if int(ch) not in o: o.append(int(ch))
        order = tuple(o)
    return order
full = '--full' in sys.argv
maxper, checked = 0, 0
for key, cyc in cycles.items():
    c0 = cyc[0]
    sup = [d for d in range(10) if c0[d]]
    if not sup or (len(sup) > 8 and not full): continue
    perms = list(permutations(sup))
    Fm = {p: turn(p, cyc) for p in perms}
    for p in perms:                                        # cycle length of the orbit's eventual cycle
        seen_p, x = {}, p
        while x not in seen_p:
            seen_p[x] = len(seen_p); x = Fm[x]
        L = len(seen_p) - seen_p[x]
        maxper = max(maxper, L * len(cyc))
    checked += 1
print('order dynamics checked on %d count cycles%s: largest number-cycle period %d'
      % (checked, '' if full else ' (supports <= 8)', maxper))
assert maxper == 4

# 4. The data files: every listed self-inventoried prime satisfies g(n) = n and is prime (all have at most 19
#    digits, below 2^64, where BPSW is exact); every listed cycle of period 2, 3 or 4 is a cycle of g;
#    exactly 4 of these cycles consist of primes only.
import os, gmpy2
here = os.path.dirname(os.path.abspath(__file__))
SP = [l.strip() for l in open(os.path.join(here, '..', 'data', 'p207_selfinv_primes.txt'))]
assert len(SP) == 9278 and all(g(x) == x and gmpy2.is_prime(int(x)) and int(x) < 2**64 for x in SP)
assert max(map(int, SP)) == 1111918171614121013
CY = [l.split()[1:] for l in open(os.path.join(here, '..', 'data', 'p207_cycles.txt'))]
assert len(CY) == 2277
for line in open(os.path.join(here, '..', 'data', 'p207_cycles.txt')):
    per, *c = line.split()
    assert int(per.rstrip(':')) == len(c) and all(g(c[i]) == c[(i + 1) % len(c)] for i in range(len(c)))
    assert len(set(c)) == len(c)
allprime = [c for c in CY if all(gmpy2.is_prime(int(x), 40) for x in c)]
assert len(allprime) == 4 and all(len(c) == 2 for c in allprime)
print('data: %d self-inventoried primes (largest %s), %d cycles of period 2-4, %d of them all prime'
      % (len(SP), max(SP, key=int), len(CY), len(allprime)))

# 5. (--full) Every self-inventoried number is a fixed point of the order map on a count fixed point:
#    766347 of them in total, 9278 prime.
if full:
    tot = 0; primes = 0
    for key, cyc in cycles.items():
        if len(cyc) != 1: continue
        c0 = cyc[0]
        sup = [d for d in range(10) if c0[d]]
        if not sup: continue
        for p in permutations(sup):
            if turn(p, cyc) == p:
                n = ''.join(str(c0[d]) + str(d) for d in p)
                assert g(n) == n
                tot += 1; primes += gmpy2.is_prime(int(n))
    assert tot == 766347 and primes == 9278
    print('all self-inventoried numbers: %d, of which %d prime' % (tot, primes))
print('all checks passed')
