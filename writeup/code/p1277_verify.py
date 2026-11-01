# Puzzle 1277 (writeup/1277.md): check the 100 x 100 grid in ../data/p1277_grid100.txt.
# It must hold every number 1..10000 exactly once; the score is the number of edge-neighbour pairs (up/down/left/right)
# that are coprime. Also prints where those pairs are, by the kinds of the two numbers, and the counts behind the
# construction: 561 'lonely' numbers (1 and the primes > 5000) are coprime to every other entry.
# Run: python3 p1277_verify.py   (a few seconds; needs sympy)
import os
from math import gcd
from collections import Counter
from sympy import isprime, primepi

here = os.path.dirname(os.path.abspath(__file__))
G = [list(map(int, l.split())) for l in open(os.path.join(here, '..', 'data', 'p1277_grid100.txt'))]
n = len(G); N = n * n
assert all(len(row) == n for row in G)
assert sorted(v for row in G for v in row) == list(range(1, N + 1))

def kind(v):
    if v == 1 or (isprime(v) and 2 * v > N): return 'lonely (1, p > 5000)'
    if isprime(v):
        for lo, name in [(N / 3, 'prime in (3333, 5000]'), (N / 4, 'prime in (2500, 3333]'), (N / 5, 'prime in (2000, 2500]'),
                         (100, 'prime in (100, 2000]')]:
            if v > lo: return name
        return 'prime < 100'
    return 'even' if v % 2 == 0 else 'odd composite'

score = 0; cat = Counter()
for r in range(n):
    for c in range(n):
        for r2, c2 in ((r, c + 1), (r + 1, c)):
            if r2 < n and c2 < n and gcd(G[r][c], G[r2][c2]) == 1:
                score += 1; cat[' & '.join(sorted([kind(G[r][c]), kind(G[r2][c2])]))] += 1
edges = 2 * n * (n - 1)
print('%d x %d grid: %d coprime neighbour pairs out of %d (%.2f %%)' % (n, n, score, edges, 100 * score / edges))
for k, v in cat.most_common(): print('   %-50s %d' % (k, v))
lonely = 1 + int(primepi(N)) - int(primepi(N // 2))
assert lonely == 561
print('lonely numbers: %d; every edge touching one is coprime' % lonely)
assert score == 1452, score
print('all checks passed')
