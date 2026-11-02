# Puzzles 764, 766, 767 (writeup/764.md): put distinct positive integers on the vertices of a graph so that the gcds of
# the edges are exactly 1, 2, ..., E (E = number of edges), with the smallest possible sum.
# Exact model (OR-Tools CP-SAT): v[x] in [1, B] all different; g[e] in [1, E] all different; for every edge (x, y):
# (v[x], v[y], g[e]) is in the table {(a, b, gcd(a, b))}. Minimise sum v[x].
# Bound: a solution with sum <= S has every value <= B = S - n(n-1)/2 (the other n-1 values are distinct and >= 1, 2, ...).
# So with S = record, the search over [1, B] decides whether the record can be matched or beaten.
# Usage: python3 gcd_graph_cpsat.py path N S [time_s]      (puzzle 764: a row of N numbers)
#        python3 gcd_graph_cpsat.py grid N S [time_s]      (puzzle 766: an N x N grid)
#        python3 gcd_graph_cpsat.py cube N S [time_s]      (puzzle 767: an N x N x N cube)
import sys, math, itertools
from ortools.sat.python import cp_model

kind, N, S = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
tl = float(sys.argv[4]) if len(sys.argv) > 4 else 600
if kind == 'path':
    cells = list(range(N)); edges = [(i, i + 1) for i in range(N - 1)]
elif kind == 'grid':
    cells = [(r, c) for r in range(N) for c in range(N)]
    edges = [((r, c), (r, c + 1)) for r in range(N) for c in range(N - 1)] + [((r, c), (r + 1, c)) for r in range(N - 1) for c in range(N)]
else:
    cells = [(a, b, c) for a in range(N) for b in range(N) for c in range(N)]
    edges = [(p, q) for p in cells for q in cells if p < q and sum(abs(p[i] - q[i]) for i in range(3)) == 1]
n, E = len(cells), len(edges)
B = S - n * (n - 1) // 2
print('%s N=%d: %d cells, %d edges, sum <= %d forces values <= %d' % (kind, N, n, E, S, B), flush=True)
table = [(a, b, math.gcd(a, b)) for a in range(1, B + 1) for b in range(1, B + 1) if math.gcd(a, b) <= E]
M = cp_model.CpModel()
v = {x: M.NewIntVar(1, B, 'v') for x in cells}
g = [M.NewIntVar(1, E, 'g') for _ in edges]
M.AddAllDifferent(list(v.values())); M.AddAllDifferent(g)
for e, (x, y) in enumerate(edges):
    M.AddAllowedAssignments([v[x], v[y], g[e]], table)
M.Add(sum(v.values()) <= S)
M.Minimize(sum(v.values()))
solver = cp_model.CpSolver(); solver.parameters.max_time_in_seconds = tl; solver.parameters.num_workers = 2
st = solver.Solve(M)
print('status', solver.StatusName(st), end=' ')
if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    vals = {x: solver.Value(v[x]) for x in cells}
    print('sum', sum(vals.values()), 'bound', solver.BestObjectiveBound())
    if kind == 'path': print(' '.join(str(vals[x]) for x in cells))
    elif kind == 'grid':
        for r in range(N): print(' '.join(str(vals[(r, c)]) for c in range(N)))
    else:
        for a in range(N):
            for b in range(N): print(' '.join(str(vals[(a, b, c)]) for c in range(N)))
            print('--')
    print('gcds', sorted(math.gcd(vals[x], vals[y]) for x, y in edges))
else:
    print('(no solution with sum <= %d in the time limit)' % S)
