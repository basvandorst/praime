# Puzzle 771 (writeup/771.md): N x N grids of digits in which every prime 2..P can be read along contiguous cells in a
# straight line (8 directions, no cell twice). Checks our grids by reading every straight segment, checks the page's
# grids the same way, and (with --sat) re-proves maximality with a SAT solver: for the next prime after P the formula
#   one digit per cell, and for each prime q <= P' an OR over all straight placements of q
# is unsatisfiable (N = 3 and 4; for N = 5 the solver does not finish). Run: python3 p771_verify.py [--sat]
# (needs sympy; --sat needs python-sat, about 2 minutes)
import sys, itertools
from sympy import primerange, nextprime

DIRS = [(0, 1), (1, 0), (1, 1), (1, -1), (0, -1), (-1, 0), (-1, -1), (-1, 1)]

def reads(G):
    N = len(G); out = set()
    for r in range(N):
        for c in range(N):
            for dr, dc in DIRS:
                s, rr, cc = '', r, c
                while 0 <= rr < N and 0 <= cc < N:
                    s += G[rr][cc]; out.add(s); rr += dr; cc += dc
    return out

def missing(G, P):
    R = reads(G)
    return [p for p in primerange(2, P + 1) if str(p) not in R]

def sat(N, P):
    from pysat.solvers import Cadical153
    top = [0]
    def new():
        top[0] += 1; return top[0]
    cells = [(r, c) for r in range(N) for c in range(N)]
    x = {(rc, d): new() for rc in cells for d in range(10)}
    cls = []
    for rc in cells:
        cls.append([x[(rc, d)] for d in range(10)])
        cls += [[-x[(rc, a)], -x[(rc, b)]] for a, b in itertools.combinations(range(10), 2)]
    for p in primerange(2, P + 1):
        w = str(p); k = len(w); opts = set()
        for (r, c) in cells:
            for dr, dc in (DIRS if k > 1 else [(0, 1)]):
                pts = tuple((r + i * dr, c + i * dc) for i in range(k))
                if all(0 <= a < N and 0 <= b < N for a, b in pts): opts.add(pts)
        ys = []
        for pts in opts:
            y = new(); ys.append(y)
            cls += [[-y, x[(rc, int(ch))]] for rc, ch in zip(pts, w)]
        cls.append(ys)
    return Cadical153(bootstrap_with=cls).solve()

PAGE = [('Kamenetsky N=4', ['6101', '4701', '3923', '1851'], 131), ('Vantieghem N=4', ['1921', '0347', '7159', '6138'], 157),
        ('Kamenetsky N=5', ['19311', '26101', '12791', '78359', '94131'], 223),
        ('Vantieghem N=5', ['10761', '93113', '24531', '17982', '19127'], 241)]
for name, G, P in PAGE:
    m = missing(G, P)
    print('%-16s claimed P = %d: %s' % (name, P, 'valid' if not m else 'INVALID, cannot read %s ...' % m[:4]))
GRIDS = {}
for line in open(__file__.replace('p771_verify.py', '../data/p771_grids.txt')):
    w = line.split()
    if w: GRIDS[int(w[0])] = (int(w[1]), w[2:])
for N, (P, G) in sorted(GRIDS.items()):
    assert len(G) == N and all(len(row) == N for row in G)
    assert not missing(G, P)
    print('ours N = %d: every prime up to %d readable: %s' % (N, P, ' / '.join(G)))
    if '--sat' in sys.argv and N <= 4:
        q = nextprime(P)
        assert not sat(N, q)
        print('   SAT: no %d x %d grid holds every prime up to %d, so %d is the maximum' % (N, N, q, P))
print('all checks passed')
