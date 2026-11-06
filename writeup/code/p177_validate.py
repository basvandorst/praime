"""Puzzle 177: 1..n^2 in an n x n square, every two horizontal or vertical neighbours sum to a prime.

1. The page's n = 4 example is checked.
2. n = 3: all 9! arrangements are tested (none works), and the reason is shown: the centre is adjacent to all four
   edge cells, which must hold 2, 4, 6, 8 (parity), and no odd number has a prime sum with each of 2, 4, 6, 8.
3. Strategy: odd numbers on the cells with i + j even, even numbers on the others (forced by parity for odd n, and
   the natural choice for even n); then a min-conflicts search: repeatedly take a cell on a composite-sum edge and
   swap it with the same-colour cell that lowers the number of composite-sum edges most (random swap when stuck);
   restart from a new random start after a move budget.  Sizes from the command line (default 4..12 and 16, 20).
"""
import sys, random, time
from itertools import permutations
from sympy import primerange, isprime

def ok(G):
    n = len(G)
    assert sorted(v for r in G for v in r) == list(range(1, n * n + 1))
    return all(isprime(G[i][j] + G[i][j + 1]) for i in range(n) for j in range(n - 1)) and \
           all(isprime(G[i][j] + G[i + 1][j]) for i in range(n - 1) for j in range(n))

# 1. page example
ex4 = [[9, 14, 5, 12], [8, 3, 2, 11], [15, 16, 1, 6], [4, 13, 10, 7]]
assert ok(ex4)
print("page's 4 x 4 example: valid")

# 2. n = 3
edges3 = [(a, b) for a in range(9) for b in range(9) if (b == a + 1 and a % 3 != 2) or b == a + 3]
count = sum(1 for P in permutations(range(1, 10)) if all(isprime(P[a] + P[b]) for a, b in edges3))
assert count == 0
partners = {e: [o for o in (1, 3, 5, 7, 9) if isprime(e + o)] for e in (2, 4, 6, 8)}
common = set.intersection(*(set(v) for v in partners.values()))
assert common == set()
print("n = 3: no arrangement among all 9! = 362880; prime-sum partners of 2, 4, 6, 8 among the odd numbers:",
      partners, "-> no common partner for the centre")

# 3. min-conflicts search
def build(n, seed, budget):
    N = n * n
    pr = [False] * (2 * N + 2)
    for p in primerange(2, 2 * N + 2):
        pr[p] = True
    rng = random.Random(seed)
    cells = [(i, j) for i in range(n) for j in range(n)]
    nb = {c: [(c[0] + di, c[1] + dj) for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1))
              if 0 <= c[0] + di < n and 0 <= c[1] + dj < n] for c in cells}
    black = [c for c in cells if (c[0] + c[1]) % 2 == 0]
    white = [c for c in cells if (c[0] + c[1]) % 2 == 1]
    odds, evens = list(range(1, N + 1, 2)), list(range(2, N + 1, 2))
    rng.shuffle(odds); rng.shuffle(evens)
    g = dict(zip(black, odds)); g.update(zip(white, evens))
    def bad(c, v):
        return sum(1 for d in nb[c] if not pr[v + g[d]])
    badc = {c: bad(c, g[c]) for c in cells}
    confl = {c for c in cells if badc[c]}
    cost = sum(badc.values()) // 2
    def refresh(c):
        for d in nb[c] + [c]:
            badc[d] = bad(d, g[d])
            if badc[d]: confl.add(d)
            else: confl.discard(d)
    moves = 0
    while cost > 0 and moves < budget:
        moves += 1
        c = rng.choice(sorted(confl))
        same = black if (c[0] + c[1]) % 2 == 0 else white
        best, bestd = None, 1
        for d in rng.sample(same, min(len(same), 150)):
            if d == c or d in nb[c]:
                continue
            vc, vd = g[c], g[d]
            before = bad(c, vc) + bad(d, vd)
            g[c], g[d] = vd, vc
            delta = bad(c, vd) + bad(d, vc) - before
            g[c], g[d] = vc, vd
            if delta < bestd:
                best, bestd = d, delta
        if best is None:
            best = rng.choice(same)
            while best == c or best in nb[c]:
                best = rng.choice(same)
        g[c], g[best] = g[best], g[c]
        refresh(c); refresh(best)
        cost = sum(badc.values()) // 2
    return ([[g[(i, j)] for j in range(n)] for i in range(n)] if cost == 0 else None), moves

sizes = [int(a) for a in sys.argv[1:]] or list(range(4, 13)) + [16, 20]
for n in sizes:
    t0 = time.time()
    for seed in range(100):
        G, moves = build(n, seed, 60 * n * n)
        if G:
            break
    assert G is not None and ok(G), n
    print(f"n = {n}: square built and checked (seed {seed}, {moves} moves, {time.time() - t0:.0f} s)", flush=True)
    if n == 5:
        for r in G: print("   ", " ".join(f"{v:3d}" for v in r))
