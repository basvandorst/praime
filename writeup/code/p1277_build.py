# Puzzle 1277 (writeup/1277.md): structured start for the 100 x 100 grid, with a mask of cells that stay fixed.
#   ring 0 (border)      : 396 'lonely' numbers (1 and primes > 5000), coprime to everything
#   ring 1               : 223 primes in (2500, 5000] with 2p directly inside (ring 2), and the other 165 lonely numbers
#   ring 2 under lonely  : remaining primes p in (2500, 3333] in groups  p 3p 3q q  (3p shares 3 with 3q), 2p inside p,
#                          a multiple of 6 inside 3p
#   interior             : primes 101..1249 as plus shapes p | 2p 4p 6p 8p;
#                          primes 1250..2500 with arms 2p, 4p, (6p or 5p or nothing) and 3p, where 3p's other three
#                          neighbours are multiples of 6 (5p's: multiples of 10)
# Everything else is placed freely; p1277_grid.c then improves the free cells by local swaps.
# Usage: python3 p1277_build.py start.txt mask.txt
import sys
from sympy import primerange, isprime

n = 100; N = n * n
G = [[0] * n for _ in range(n)]; F = [[0] * n for _ in range(n)]
used = set()
def put(r, c, v, frozen=1):
    assert G[r][c] == 0 and v not in used and 1 <= v <= N, (r, c, v)
    G[r][c] = v; F[r][c] = frozen; used.add(v)

P = list(primerange(2, N + 1))
lonely = [1] + [p for p in P if 2 * p > N]; lonely.sort(reverse=True)
L1 = sorted([p for p in P if 3 * p > N and 2 * p <= N], reverse=True)          # one multiple: 2p
L2 = sorted([p for p in P if 4 * p > N and 3 * p <= N], reverse=True)          # 2p, 3p
mid = sorted([p for p in P if 1250 <= p <= 2500])                               # 2p, 4p, 3p (+ 6p or 5p)
small = [p for p in P if 101 <= p < 1250]                                        # 2p, 4p, 6p, 8p
pool6 = [v for v in range(6, N + 1, 6)]                                          # multiples of 6, taken as needed
pool10 = [v for v in range(10, N + 1, 10) if v % 6]                              # multiples of 10 (not of 6)
reserved = set()
for p in L1 + L2: reserved |= {2 * p, 3 * p}
for p in mid: reserved |= {k * p for k in (2, 3, 4, 5, 6) if k * p <= N}
for p in small: reserved |= {2 * p, 4 * p, 6 * p, 8 * p}
def take(pool, avoid):
    while pool[-1] in used or pool[-1] in avoid or pool[-1] in reserved: pool.pop()
    return pool.pop()

def ring_cells(k):
    lo, hi = k, n - 1 - k; out = []
    for c in range(lo, hi): out.append((lo, c))
    for r in range(lo, hi): out.append((r, hi))
    for c in range(hi, lo, -1): out.append((hi, c))
    for r in range(hi, lo, -1): out.append((r, lo))
    return out
def inward(r, c, k):
    lo, hi = k, n - 1 - k
    if r == lo and lo < c < hi: return (r + 1, c)
    if r == hi and lo < c < hi: return (r - 1, c)
    if c == lo and lo < r < hi: return (r, c + 1)
    if c == hi and lo < r < hi: return (r, c - 1)
    return None                                                                   # a corner

# ring 0: lonely
it = iter(lonely)
for r, c in ring_cells(0): put(r, c, next(it))
rest_lonely = list(it)
# ring 1: primes with 2p inside first, lonely at the end (and on the corners)
semi = L1 + L2[:388 - len(rest_lonely)]
cells1 = ring_cells(1)
corner2 = {(2, 2), (2, n - 3), (n - 3, 2), (n - 3, n - 3)}
cand = [rc for rc in cells1 if inward(*rc, 1) is not None and inward(*rc, 1) not in corner2]
nsemi = 388 - len(rest_lonely)
for (r, c), p in zip(cand[:nsemi], semi[:nsemi]):
    put(r, c, p); put(*inward(r, c, 1), 2 * p)
lon_cells = [rc for rc in cells1 if G[rc[0]][rc[1]] == 0]
assert len(lon_cells) == len(rest_lonely)
for (r, c), v in zip(lon_cells, rest_lonely): put(r, c, v)
# ring 2 cells under the lonely arc: groups p 3p 3q q of the remaining L2 primes
left2 = [p for p in L2 if p not in used]
arc = [inward(r, c, 1) for (r, c) in cand[nsemi:]]
arc = [rc for rc in arc if rc is not None and G[rc[0]][rc[1]] == 0]
i = 0
while i + 4 <= len(arc) and len(left2) >= 2:
    p, q = left2.pop(), left2.pop()
    for (r, c), v in zip(arc[i:i + 4], [p, 3 * p, 3 * q, q]):
        put(r, c, v)
        inn = inward(r, c, 2)
        if inn is not None and G[inn[0]][inn[1]] == 0: put(*inn, 2 * v if v in (p, q) else take(pool6, set()))
    i += 4
leftover_L2 = left2
# interior structures, rows/cols 4..95, in boxes
boxes = []
r0 = 4
def plus_box(r, c, p):            # 3x3 box with top-left (r, c): p at centre, 2p 4p 6p 8p around
    put(r + 1, c + 1, p); put(r, c + 1, 2 * p); put(r + 1, c, 4 * p); put(r + 1, c + 2, 6 * p); put(r + 2, c + 1, 8 * p)
def mid_box(r, c, p):             # 4x3 box: p centre (r+1, c+1); 2p up, 4p left, right 6p/5p/free, 3p down with 6m around
    put(r + 1, c + 1, p); put(r, c + 1, 2 * p); put(r + 1, c, 4 * p)
    if 6 * p <= N: put(r + 1, c + 2, 6 * p)
    elif 5 * p <= N:
        put(r + 1, c + 2, 5 * p)
    put(r + 2, c + 1, 3 * p)
    for rc in ((r + 3, c + 1), (r + 2, c), (r + 2, c + 2)): put(*rc, take(pool6, {6 * p}))
row, col = 4, 4
queue = [('mid', p) for p in mid] + [('plus', p) for p in small]
# mid boxes are 4 rows x 3 cols; plus boxes 3 x 3 (packed with 1 free column between boxes)
for kind, p in queue:
    h = 4 if kind == 'mid' else 3
    if col + 3 > 96: col = 4; row += 5 if kind == 'mid' else 4
    if row + h > 96: raise SystemExit('interior full')
    (mid_box if kind == 'mid' else plus_box)(row, col, p)
    col += 4
# 5p cells need multiples of 10 around them: mark which; leave the rest to the annealing (it sees 5p fixed)
# everything else: free
rest = [v for v in range(1, N + 1) if v not in used]
free = [(r, c) for r in range(n) for c in range(n) if G[r][c] == 0]
assert len(rest) == len(free)
rest.sort(key=lambda v: (v % 2, v))            # evens first (they fill the space between structures)
for (r, c), v in zip(free, rest): G[r][c] = v; F[r][c] = 0
open(sys.argv[1], 'w').write('\n'.join(' '.join(map(str, row_)) for row_ in G) + '\n')
open(sys.argv[2], 'w').write('\n'.join(' '.join(map(str, row_)) for row_ in F) + '\n')
print('frozen cells', sum(map(sum, F)), ' leftover primes in (2500, 3333]:', leftover_L2, ' last box row', row)
