# Puzzle 685 Q3: Hamiltonian paths of 1..N with every adjacent sum odd and composite, for all N >= 14.
# Construction: blocks 30m + r arranged by a fixed residue order; every adjacent sum is odd and divisible by 3 or 5.
import sys
sys.setrecursionlimit(10000)
ok_res = lambda s: s % 2 == 1 and (s % 3 == 0 or s % 5 == 0)       # residue condition (mod 30 is enough)
def dfs_path(nodes, start_ok, end_ok, edge_ok):
    nodes = list(nodes); n = len(nodes)
    order = []
    used = set()
    def go():
        if len(order) == n: return end_ok(order[-1])
        last = order[-1]
        for v in sorted(nodes, key=lambda v: (v * 7) % 31):
            if v not in used and edge_ok(last, v):
                used.add(v); order.append(v)
                if go(): return True
                used.discard(v); order.pop()
        return False
    for s in nodes:
        if start_ok(s):
            used.add(s); order.append(s)
            if go(): return list(order)
            used.discard(s); order.pop()
    return None
# 1. the periodic block: an order of residues 1..30 with every internal sum and the wrap r_30 + r_1 allowed
edge = lambda a, b: ok_res(a + b)
block = None
for first in range(1, 31):
    block = dfs_path(range(1, 31), lambda s: s == first, lambda last: ok_res(last + first), edge)
    if block: break
print("block:", block)
# 2. a tail for every rho = N mod 30: order of residues {1..rho} starting with a v compatible with r_30 (v + r_30 allowed)
tails = {0: []}
for rho in range(1, 30):
    t = dfs_path(range(1, rho + 1), lambda s: ok_res(s + block[-1]), lambda last: True, edge)
    tails[rho] = t
    print("rho", rho, "tail", t)
missing = [r for r, t in tails.items() if t is None]
print("rho without a tail:", missing)

# 3. long tails for rho = 4..12: offsets {1..30+rho} after a full block, start compatible with block[-1]
for rho in missing:
    t = dfs_path(range(1, 31 + rho), lambda s: ok_res(s + block[-1]), lambda last: True, edge)
    tails[rho] = ('long', t); print("rho", rho, "long tail", t)
assert all(tails[r][1] is not None for r in missing)
# 4. base cases 14 <= N < 120 by direct search with real sums (odd, composite)
from sympy import isprime
comp_odd = lambda s: s % 2 == 1 and s > 1 and not isprime(s)
base = {}
for N in range(14, 120):
    p = dfs_path(range(1, N + 1), lambda s: True, lambda last: True, lambda a, b: comp_odd(a + b))
    assert p, N; base[N] = p
print("base cases 14..119 found")
def build(N):
    if N < 120: return base[N]
    M, rho = divmod(N, 30)
    path = []
    if isinstance(tails[rho], tuple):          # long tail after blocks 0..M-2
        for m in range(M - 1): path += [30 * m + r for r in block]
        path += [30 * (M - 1) + r for r in tails[rho][1]]
    else:
        for m in range(M): path += [30 * m + r for r in block]
        path += [30 * M + r for r in tails[rho]]
    return path
for N in range(14, 3001):
    p = build(N)
    assert sorted(p) == list(range(1, N + 1)), N
    assert all(comp_odd(a + b) for a, b in zip(p, p[1:])), N
print("construction verified for every N in 14..3000")
import json
json.dump({"block": block, "tails": {str(k): v for k, v in tails.items()}}, open("construction685.json", "w"))
