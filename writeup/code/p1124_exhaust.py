"""Puzzles 1123 / 1124, Python version of p1124_exhaust.c.  Triangle with R rows holding 0..R(R+1)/2-1; along
every straight line (rows and both diagonal directions) every contiguous segment is read forwards and backwards,
the numbers concatenated in decimal; count the distinct primes.  Exhaustive DFS with the symmetry broken by
requiring top corner < left corner < right corner.  Feasible in Python for R = 4 (tetractys: maximum 24, nine
arrangements); R = 5 needs the C program.      Usage: p1124_exhaust.py R T"""
import sys
from gmpy2 import is_prime
R, T = int(sys.argv[1]), int(sys.argv[2])
NC = R * (R + 1) // 2
cell = {}; k = 0
for r in range(R):
    for i in range(r + 1): cell[(r, i)] = k; k += 1
lines = [[cell[(r, i)] for i in range(r + 1)] for r in range(R)]
lines += [[cell[(r, i)] for r in range(i, R)] for i in range(R)]
lines += [[cell[(r, r - j)] for r in range(j, R)] for j in range(R)]
segs = []
for ln in lines:
    for a in range(len(ln)):
        for b in range(a, len(ln)):
            segs.append(tuple(ln[a:b + 1]))
            if b > a: segs.append(tuple(ln[b:a - 1 if a else None:-1]) if a else tuple(ln[b::-1]))
segs = list(set(segs))
order = [cell[(0, 0)], cell[(R - 1, 0)], cell[(R - 1, R - 1)]] + [c for c in range(NC) if c not in (cell[(0, 0)], cell[(R - 1, 0)], cell[(R - 1, R - 1)])]
depth = {c: d for d, c in enumerate(order)}
complete_at = [[] for _ in range(NC)]
for s in segs: complete_at[max(depth[c] for c in s)].append(s)
memo = {}
def prime_value(vals):
    if vals in memo: return memo[vals]
    v = int("".join(map(str, vals))); r = v if is_prime(v) else None
    memo[vals] = r; return r
val = [-1] * NC; used = [False] * NC; cnt = {}; results = []; best = [0]
def dfs(d, distinct):
    if d == NC:
        best[0] = max(best[0], distinct)
        if distinct >= T: results.append((distinct, val[:]))
        return
    c = order[d]
    for v in range(NC):
        if used[v]: continue
        if c == order[1] and v < val[order[0]]: continue
        if c == order[2] and (v < val[order[0]] or v < val[order[1]]): continue
        val[c] = v; used[v] = True; added = []
        for s in complete_at[d]:
            p = prime_value(tuple(val[x] for x in s))
            if p is not None:
                cnt[p] = cnt.get(p, 0) + 1
                if cnt[p] == 1: added.append(p)
                else: added.append(None)
        dfs(d + 1, distinct + sum(1 for a in added if a is not None))
        for s in complete_at[d]:
            p = prime_value(tuple(val[x] for x in s))
            if p is not None: cnt[p] -= 1
        used[v] = False; val[c] = -1
dfs(0, 0)
for d, arr in results: print(d, "primes:", " ".join(map(str, arr)))
print(f"R={R}: maximum {best[0]}, arrangements with >= {T} primes: {len(results)}")
