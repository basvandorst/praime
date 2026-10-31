# Puzzle 120, Q4: no n x n "full torus" primal square exists for n = 4, 6, 8 (writeup/120.md).
# A full-torus primal square is an arrangement of 1..n^2 where every two grid neighbours and every two mirror cells
# (r,c)-(r,n-1-c) and (r,c)-(n-1-r,c) sum to a prime (Jobling's definition on the page).
# Run: python3 p120_verify.py          (n = 4 and 6, about 10 seconds)
#      python3 p120_verify.py --full   (also n = 8, a few minutes; the C version p120_mod3.c does it in seconds)
import sys
from itertools import permutations
from sympy import isprime


def graph(n, rivera=False):
    """Edges of the constraint graph. rivera=True drops mirror pairs between main-diagonal cells (Rivera's wording)."""
    E = set()
    for r in range(n):
        for c in range(n):
            a = (r, c)
            if c + 1 < n: E.add(frozenset([a, (r, c + 1)]))
            if r + 1 < n: E.add(frozenset([a, (r + 1, c)]))
            diag = (r == c or r + c == n - 1)
            if rivera and diag:
                continue
            if c != n - 1 - c: E.add(frozenset([a, (r, n - 1 - c)]))
            if r != n - 1 - r: E.add(frozenset([a, (n - 1 - r, c)]))
    return [tuple(e) for e in E]


def colourable(n, orient, rivera=False):
    """Is there an assignment of residues mod 3 to the cells that a primal square would induce?
    Rules: neighbours have opposite parity (cells with (r+c) % 2 == 0 hold numbers of parity `orient`);
    no neighbours 0-0 and no neighbours 1-2, except one edge between the cells of 1 (odd, residue 1) and
    2 (even, residue 2); the number of cells of each (parity, residue) equals that of the numbers 1..n^2."""
    cells = [(r, c) for r in range(n) for c in range(n)]           # row-major order (the C program uses another)
    nb = {x: [] for x in cells}
    for a, b in graph(n, rivera):
        nb[a].append(b); nb[b].append(a)
    par = {(r, c): (orient if (r + c) % 2 == 0 else 1 - orient) for r, c in cells}   # 0 = even number, 1 = odd
    need = {(p, x): 0 for p in (0, 1) for x in (0, 1, 2)}
    for v in range(1, n * n + 1):
        need[(v % 2, v % 3)] += 1
    have = dict.fromkeys(need, 0)
    col = {}
    state = {'ex': 0, 'nodes': 0}

    def bt(k):
        state['nodes'] += 1
        if k == len(cells):
            return True
        a = cells[k]
        p = par[a]
        for x in (0, 1, 2):
            if have[(p, x)] >= need[(p, x)]:
                continue
            ex = 0; ok = True
            for b in nb[a]:
                if b not in col:
                    continue
                y = col[b]
                if x == 0 and y == 0:
                    ok = False; break
                if {x, y} == {1, 2}:
                    one_cell_odd = (p == 1) if x == 1 else (par[b] == 1)
                    two_cell_even = (p == 0) if x == 2 else (par[b] == 0)
                    if one_cell_odd and two_cell_even:
                        ex += 1
                    else:
                        ok = False; break
            if not ok or state['ex'] + ex > 1:
                continue
            col[a] = x; have[(p, x)] += 1; state['ex'] += ex
            if bt(k + 1):
                return True
            del col[a]; have[(p, x)] -= 1; state['ex'] -= ex
        return False

    return bt(0), state['nodes']


# 1. The two parity orientations, Jobling's definition: no residue colouring exists, so no square exists.
sizes = [4, 6] + ([8] if '--full' in sys.argv else [])
for n in sizes:
    for orient in (0, 1):
        ok, nodes = colourable(n, orient)
        print('n = %d, orientation %d: residue colouring %s (%d nodes)' % (n, orient, 'exists' if ok else 'impossible', nodes))
        assert not ok

# 2. Rivera's weaker wording (mirror pairs between diagonal cells dropped): still impossible for n = 4 and 6.
for n in (4, 6):
    for orient in (0, 1):
        assert not colourable(n, orient, rivera=True)[0]

# 3. n = 4 by hand: 14 has only the partners 3, 5, 9, 15 (14 + x prime), and for the pairs {5, x} the only even
#    common partners are 2 and 8 - three 4-cycles through 14 would need three distinct such cells.
partners = {e: [o for o in range(1, 17, 2) if isprime(e + o)] for e in range(2, 17, 2)}
assert partners[14] == [3, 5, 9, 15]
for x in (3, 9, 15):
    assert [e for e in range(2, 17, 2) if e != 14 and isprime(e + 5) and isprime(e + x)] == [2, 8]

# 4. Direct exhaustive check for n = 4 (both readings): no arrangement of 1..16 at all.
def count4(rivera):
    E = graph(4, rivera)
    nb = {}
    for a, b in E:
        nb.setdefault(a, []).append(b); nb.setdefault(b, []).append(a)
    cells = [(r, c) for r in range(4) for c in range(4)]
    val = {}
    used = set()
    sols = 0
    def bt(k):
        nonlocal sols
        if k == 16:
            sols += 1; return
        a = cells[k]
        for x in range(1, 17):
            if x in used: continue
            if all(isprime(x + val[b]) for b in nb[a] if b in val):
                val[a] = x; used.add(x); bt(k + 1); used.discard(x); del val[a]
    bt(0)
    return sols
assert count4(False) == 0 and count4(True) == 0
print('all checks passed')
