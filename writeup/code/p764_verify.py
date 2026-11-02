# Puzzles 764, 766, 767 (writeup/764.md): distinct positive integers on a row / grid / cube whose neighbour gcds are
# exactly 1..E, with minimal sum. Checks every arrangement on the page and ours (distinct values, gcd set = 1..E, sum),
# and the bound used by the exact searches: with sum <= S and n distinct values, every value is <= S - n(n-1)/2.
# Run: python3 p764_verify.py   (a second; no dependencies)
import os
from math import gcd

def edges(kind, N):
    if kind == 'path': return [(i, i + 1) for i in range(N - 1)], N
    if kind == 'grid':
        idx = lambda r, c: r * N + c
        return [(idx(r, c), idx(r, c + 1)) for r in range(N) for c in range(N - 1)] + \
               [(idx(r, c), idx(r + 1, c)) for r in range(N - 1) for c in range(N)], N * N
    idx = lambda a, b, c: (a * N + b) * N + c
    E = []
    for a in range(N):
        for b in range(N):
            for c in range(N):
                if a + 1 < N: E.append((idx(a, b, c), idx(a + 1, b, c)))
                if b + 1 < N: E.append((idx(a, b, c), idx(a, b + 1, c)))
                if c + 1 < N: E.append((idx(a, b, c), idx(a, b, c + 1)))
    return E, N ** 3

def check(name, kind, N, vals):
    E, n = edges(kind, N)
    assert len(vals) == n and len(set(vals)) == n and min(vals) >= 1, name
    gs = sorted(gcd(vals[x], vals[y]) for x, y in E)
    assert gs == list(range(1, len(E) + 1)), (name, gs)
    S = sum(vals)
    print('%-34s sum %7d  (%d values, gcds 1..%d; any solution with this sum has values <= %d)' % (name, S, n, len(E), S - n * (n - 1) // 2))
    return S

P = {}
P['764 N=10 (page, 102)'] = ('path', 10, [3, 9, 18, 12, 8, 16, 5, 10, 14, 7])
P['764 N=20 (page, 508)'] = ('path', 20, [11, 77, 14, 28, 32, 16, 40, 10, 15, 30, 12, 36, 18, 9, 39, 13, 19, 38, 34, 17])
P['764 N=50 (page, 4036)'] = ('path', 50, [int(x) for x in '43 86 94 47 41 123 93 124 116 203 49 98 28 84 42 105 70 50 75 45 90 60 40 80 32 96 48 24 36 72 54 27 117 78 26 143 33 66 44 88 136 34 323 38 114 138 46 115 185 37'.split()])
P['764 N=100 (page, 22062)'] = ('path', 100, [int(x) for x in '79 316 268 469 427 671 583 424 184 92 138 69 989 86 516 348 87 58 406 196 98 245 70 210 60 180 90 135 81 162 108 72 144 96 192 64 160 80 120 168 56 252 84 126 63 231 154 110 165 99 198 132 88 176 304 76 114 57 95 380 100 150 75 325 130 52 156 78 273 91 481 74 666 558 62 93 465 255 204 68 170 85 799 94 282 426 355 415 166 178 89 97 291 177 590 410 82 369 657 73'.split()])
P['766 3x3 (page, 251)'] = ('grid', 3, [10, 44, 11, 40, 24, 36, 35, 42, 9])
P['766 4x4 (page, 1647)'] = ('grid', 4, [187, 88, 110, 46, 102, 48, 80, 115, 117, 72, 60, 76, 182, 126, 105, 133])
P['766 5x5 (page, 6480)'] = ('grid', 5, [185, 175, 150, 460, 437, 518, 105, 120, 160, 38, 28, 252, 216, 352, 304, 68, 306, 297, 286, 377, 527, 186, 429, 156, 348])
P['767 2x2x2 (page, 338)'] = ('cube', 2, [66, 12, 63, 36, 55, 10, 56, 40])
for name, (kind, N, vals) in P.items(): check(name, kind, N, vals)
here = os.path.dirname(os.path.abspath(__file__))
for line in open(os.path.join(here, '..', 'data', 'p764_solutions.txt')):
    w = line.split()
    if not w or w[0].startswith('#'): continue
    kind, N, claim = w[0], int(w[1]), int(w[2]); vals = [int(x) for x in w[3:]]
    S = check('ours %s N=%d' % (kind, N), kind, N, vals)
    assert S == claim
print('all checks passed')
