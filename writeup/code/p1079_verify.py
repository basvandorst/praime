"""Puzzles 1079/1085 (prime stepping stones): independent checks.

1. Replay: in a configuration, every prime must equal the sum of its king-neighbours that have a smaller value
   (the 1, the stones of value 2 and the smaller primes), and the primes must be 3, 5, ..., P without gaps.
2. Exhaustive search for n <= NMAX stones, written independently of p1079_stones.c: primes are placed in increasing
   order; a cell for q needs neighbour sum S0 <= q with q - S0 = 2k, k new stones put on empty neighbours of the cell
   that touch no earlier prime (a stone placed in advance next to an earlier prime would have changed that prime's sum).
Usage: p1079_verify.py NMAX   (NMAX = 4 takes about 20 s, 5 about 5 minutes)
"""
import sys
from itertools import combinations
from sympy import primerange

def replay(grid):
    cells = {(r, c): int(v) for r, row in enumerate(grid) for c, v in enumerate(row.split()) if v != '.'}
    ps = sorted(v for v in cells.values() if v > 2)
    assert ps == list(primerange(3, ps[-1] + 1)), "primes not consecutive"
    assert sum(1 for v in cells.values() if v == 1) == 1
    for (r, c), v in cells.items():
        if v > 2:
            s = sum(cells.get((r + dr, c + dc), 0) for dr in (-1, 0, 1) for dc in (-1, 0, 1)
                    if (dr or dc) and cells.get((r + dr, c + dc), 0) < v)
            assert s == v, ((r, c), v, s)
    return ps[-1], sum(1 for v in cells.values() if v == 2)

NB = [(dr, dc) for dr in (-1, 0, 1) for dc in (-1, 0, 1) if dr or dc]
def search(n, P):
    best = [0]
    def rec(i, val, primecells, left):
        if i > 0: best[0] = max(best[0], P[i - 1])
        if i == len(P): return
        q = P[i]
        frontier = {(r + dr, c + dc) for (r, c) in val for dr, dc in NB} - set(val)
        for (r, c) in frontier:
            if i == 0 and (r, c) not in ((0, 1), (1, 1)): continue      # symmetry of the first move
            nbv = [val.get((r + dr, c + dc), 0) for dr, dc in NB]
            if not any(v % 2 for v in nbv): continue
            s0 = sum(nbv)
            if s0 > q or (q - s0) % 2: continue
            k = (q - s0) // 2
            if k > left: continue
            free = [(r + dr, c + dc) for dr, dc in NB if (r + dr, c + dc) not in val
                    and not any((r + dr + a, c + dc + b) in primecells for a, b in NB)]
            for S in combinations(free, k):
                v2 = dict(val); v2[(r, c)] = q
                for x in S: v2[x] = 2
                rec(i + 1, v2, primecells | {(r, c)}, left - k)
    rec(0, {(0, 0): 1}, frozenset(), n)
    return best[0]

grids = {   # ours (from p1079_stones.c) and the page's records
    4: [". 31 . . . .", ". 2 29 23 17 .", "5 41 2 2 2 13", ". 3 . 7 11 .", ". . 1 19 37 ."],
    6: [". 23 . . . . . .", "19 2 2 . . . . .", ". 17 47 11 . . . .", ". . 13 2 7 . . 61",
        ". . . 29 5 53 59 2", ". . 31 2 37 1 3 .", ". . . . . 41 . .", ". . . 2 43 . . ."],
    7: [". . 19 2 . . .", ". . 53 17 73 . .", ". . 2 13 41 . .", "2 2 37 11 . . .", "47 43 2 7 23 . 31",
        ". . 59 5 . 29 2", ". 61 2 67 1 3 .", ". . . . 71 . ."],
    8: [". . 73 2 19 . . . .", ". 71 . . 17 . . . .", "2 67 2 23 2 13 29 . .", ". 61 2 31 . 59 11 5 7",
        ". 53 2 37 . 1 3 2 .", ". 47 2 41 79 83 89 . .", ". . 43 . . . . . ."],
}
for n, g in grids.items():
    P, st = replay(g)
    assert st == n
    print(f"configuration with {n} stones reaches {P}: valid")
import os, re
here = os.path.dirname(os.path.abspath(__file__))
L = [l.strip() for l in open(os.path.join(here, '..', '..', 'data', 'puzzles', 'puzzle-1079.txt')) if l.strip()]
for i, l in enumerate(L):                       # the page's configurations (Kamenetsky), parsed from the page text
    m = re.match(r'^n=(\d+) score (\d+)$', l)
    if not m: continue
    g = []; j = i + 1
    while j < len(L) and re.match(r'^[\d.]+( [\d.]+)+$', L[j]): g.append(L[j]); j += 1
    P, st = replay(g)
    assert (st, P) == (int(m.group(1)), int(m.group(2)))
    print(f"page configuration n = {st}: reaches {P}, valid")
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 4
P = list(primerange(3, 200))
expect = {1: 13, 2: 23, 3: 31, 4: 41, 5: 53}
for n in range(1, NMAX + 1):
    b = search(n, P)
    print(f"n = {n}: largest prime reachable = {b} (independent search)")
    assert b == expect[n]
print("all checks passed")
