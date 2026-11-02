# Puzzle 764 (writeup/764.md): second, independent program for "no row of N cells with sum <= S".
# Order-first search: enumerate orders of the labels 1..N-1 by depth-first search, pruned by the exact relaxation
#   partial cost + rem(unused, last) > S, where partial cost = g1 + lcm(g1,g2) + ... + lcm(g_{k-1},g_k)
# (cells counted at their lcm), and for every complete order solve the value problem exactly: cell i = lcm_i * t_i
# with t_i >= 1, values distinct, gcd of neighbours exact, sum <= S (a small search, since lcm_i * t_i <= B).
# Usage: python3 p764_orders.py N S
import sys
from math import gcd, lcm
from functools import lru_cache

N, S = int(sys.argv[1]), int(sys.argv[2]); m = N - 1; B = S - N * (N - 1) // 2
full = (1 << m) - 1

@lru_cache(maxsize=None)
def rem(mask, g):                       # least lcm-cost of the cells after a cell whose right label is g
    if mask == 0: return 0
    return min(lcm(g, h) + rem(mask & ~(1 << (h - 1)), h) for h in range(1, m + 1) if mask >> (h - 1) & 1)

def values(order):                      # exact: can the label order carry distinct values with exact gcds, sum <= S?
    L = [order[0]] + [lcm(order[i], order[i + 1]) for i in range(m - 1)] + [order[-1]]
    used = set()
    def go(i, prev, sum_):
        if i == N: return []
        base = L[i]
        rest = sum(L[j] for j in range(i + 1, N))        # each later cell at least its lcm
        w = base
        while sum_ + w + rest <= S and w <= B:
            if w not in used and (i == 0 or gcd(prev, w) == order[i - 1]):
                used.add(w)
                r = go(i + 1, w, sum_ + w)
                used.discard(w)
                if r is not None: return [w] + r
            w += base
        return None
    return go(0, 0, 0)

orders = [0]
def dfs(order, mask, cost):
    last = order[-1]
    if mask == 0:
        orders[0] += 1
        r = values(order)
        return r
    for h in range(1, m + 1):
        if mask >> (h - 1) & 1:
            c = cost + lcm(last, h)
            nm = mask & ~(1 << (h - 1))
            if c + rem(nm, h) <= S:
                r = dfs(order + [h], nm, c)
                if r: return r
    return None

for g in range(1, m + 1):
    nm = full & ~(1 << (g - 1))
    if g + rem(nm, g) <= S:
        r = dfs([g], nm, g)
        if r:
            print('N = %d, sum <= %d: row %s sum %d (%d complete orders tried)' % (N, S, r, sum(r), orders[0])); sys.exit()
print('N = %d, sum <= %d: none (%d complete label orders examined)' % (N, S, orders[0]))
