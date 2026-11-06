"""Puzzle 176: 1..n on a circle, every two neighbours and every two opposite numbers sum to a prime (n = 2m, m odd).

1. Checks the page's n = 10 example and Jud McCranie's n = 202 arrangement.
2. The zigzag 1, n, 3, n-2, 5, ..., n-1, 2 (odd position j holds j, even position j holds n+2-j) has neighbour sums
   n+1, n+3 (and 2+1 = 3) and opposite sums m+2, 3m+2 only; the rotated zigzag (even position j holds 4-j mod n)
   has sums 3, 5, n+3, n+5, m+4, 3m+4.  Both are verified, and all n <= 3000 they cover are listed.
3. General strategy: fill position i and its opposite i+m together as an oriented pair (x_i, y_i) with x_i + y_i,
   x_i + x_{i-1} and y_i + y_{i-1} prime, by depth-first search with Warnsdorff ordering (fewest onward options
   first); the last pair must also close both neighbour sums with the first pair.  A seed gets 3000 nodes, then
   the search restarts with the next seed.  Run for every
   n = 2 mod 4 up to NMAX (argument, default 402) and for the extra sizes listed; every circle is checked.
"""
import re, sys, time, random
from sympy import isprime, nextprime, primerange

def check(circle):
    n = len(circle); m = n // 2
    if sorted(circle) != list(range(1, n + 1)):
        return False
    return all(isprime(circle[i] + circle[(i + 1) % n]) for i in range(n)) and \
           all(isprime(circle[i] + circle[i + m]) for i in range(m))

# 1. the page's examples ------------------------------------------------------
ex10 = [4, 9, 2, 1, 10, 3, 8, 5, 6, 7]                 # the drawn circle, read around
assert check(ex10)
txt = open('data/puzzles/puzzle-0176.txt').read()
seg = txt[txt.rindex('solution for n=202'):]
seg = seg[:seg.index('***')]
jud = [int(t) for t in re.findall(r'\d+', seg)][1:]
assert len(jud) == 202 and check(jud)
print("page n=10 example and Jud McCranie's n=202 circle: valid")

# 2. the zigzag families ---------------------------------------------------------
def zigzag(n, c):
    a = [0] * n
    for j in range(1, n + 1):
        a[j - 1] = j if j % 2 else ((c - j) % n or n)
    return a

def sums(a):
    n = len(a); m = n // 2
    return sorted({a[i] + a[(i + 1) % n] for i in range(n)} | {a[i] + a[i + m] for i in range(m)})

fam1, fam2 = [], []
for n in range(10, 3001, 4):          # n = 6 is degenerate (3m+4 does not occur)
    m = n // 2
    assert sums(zigzag(n, n + 2)) == sorted({3, n + 1, n + 3, m + 2, 3 * m + 2})
    assert sums(zigzag(n, 4)) == sorted({3, 5, n + 3, n + 5, m + 4, 3 * m + 4})
    if all(isprime(x) for x in (n + 1, n + 3, m + 2, 3 * m + 2)):
        assert check(zigzag(n, n + 2)); fam1.append(n)
    if all(isprime(x) for x in (n + 3, n + 5, m + 4, 3 * m + 4)):
        assert check(zigzag(n, 4)); fam2.append(n)
assert zigzag(10, 12) == [1, 10, 3, 8, 5, 6, 7, 4, 9, 2] and sums(ex10) == [3, 7, 11, 13, 17]
print("zigzag n+2 (sums 3, n+1, n+3, m+2, 3m+2), n <= 3000:", fam1)
print("zigzag 4   (sums 3, 5, n+3, n+5, m+4, 3m+4), n <= 3000:", fam2)

# 3. general strategy ----------------------------------------------------------------
def build(n, seed=0, budget=3000):
    """Positions 1..m and their opposites m+1..2m are filled together as oriented pairs (x_i, y_i):
    x_i + x_{i-1}, y_i + y_{i-1} and x_i + y_i prime; at the end x_m + y_1 and y_m + x_1 prime.
    Depth-first search, candidates ordered by fewest onward options (Warnsdorff), ties broken by the seed."""
    m = n // 2
    pr = [False] * (2 * n + 2)
    for p in primerange(2, 2 * n + 2):
        pr[p] = True
    partners = [[u for u in range(1, n + 1) if u != v and pr[v + u]] for v in range(n + 1)]
    rng = random.Random(seed)
    used = [False] * (n + 1)
    xs, ys = [1], [None]
    used[1] = True
    nodes = 0
    def onward(v):                       # number of unused u with v + u prime
        return sum(1 for u in partners[v] if not used[u])
    def dfs():
        nonlocal nodes
        nodes += 1
        if nodes > budget:
            return False
        i = len(xs)                      # next index to fill is i (0-based), current last is i-1
        x0, y0 = xs[-1], ys[-1]
        if y0 is None:                   # choose y_1 (opposite of 1)
            cand = [w for w in partners[1] if not used[w]]
            rng.shuffle(cand)
            for w in cand:
                used[w] = True; ys[0] = w
                if dfs():
                    return True
                used[w] = False; ys[0] = None
            return False
        if i == m:
            return pr[x0 + ys[0]] and pr[y0 + xs[0]]
        cand = [(v, w) for v in partners[x0] if not used[v]
                for w in partners[y0] if not used[w] and w != v and pr[v + w]]
        rng.shuffle(cand)
        scored = []
        for v, w in cand:
            used[v] = used[w] = True
            scored.append((onward(v) + onward(w), v, w))
            used[v] = used[w] = False
        scored.sort()
        for _, v, w in scored:
            used[v] = used[w] = True
            xs.append(v); ys.append(w)
            if dfs():
                return True
            xs.pop(); ys.pop()
            used[v] = used[w] = False
        return False
    sys.setrecursionlimit(max(10000, 4 * m))
    if dfs():
        return xs + ys, nodes
    return None, nodes

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 402
EXTRA = [int(a) for a in sys.argv[2:]] or [750, 1002, 2002]
hard = []
t0 = time.time()
for n in list(range(10, NMAX + 1, 4)) + EXTRA:      # n = 6 has no solution (page, Q3)
    circle = None
    for seed in range(300):
        circle, nodes = build(n, seed)
        if circle:
            break
    assert circle is not None and check(circle), n
    if seed > 0 or nodes > n:
        hard.append((n, seed, nodes))
    if n in EXTRA or n == NMAX:
        print(f"n = {n}: circle found and checked (seed {seed}, {nodes} search nodes, {time.time() - t0:.0f} s total)")
        if n == 202:
            print("   n=202:", circle)
print(f"all n = 2 (mod 4) with 10 <= n <= {NMAX} and n in {EXTRA}: circles built and checked")
print("cases that needed backtracking (n, seed, nodes):", hard[:20], "count", len(hard))
