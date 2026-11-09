# Puzzle 772, Q1: N x N grid A of distinct non-negative integers, B = sums over king neighbourhoods (cell included),
# all B distinct primes, minimise sum(B). Usage: p772_q1_cpsat.py N U TIMELIMIT   (N = 3: OPTIMAL 291 in about 10 s)
import sys
from ortools.sat.python import cp_model
from sympy import primerange
N, U, TL = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
cells = [(r, c) for r in range(N) for c in range(N)]
nb = {x: [(r, c) for r in range(x[0]-1, x[0]+2) for c in range(x[1]-1, x[1]+2) if 0 <= r < N and 0 <= c < N] for x in cells}
w = {x: sum(1 for y in cells if x in nb[y]) for x in cells}
m = cp_model.CpModel()
A = {x: m.NewIntVar(0, U, f"a{x}") for x in cells}
m.AddAllDifferent(list(A.values()))
P = list(primerange(3, 9 * U + 10))
B = {x: m.NewIntVarFromDomain(cp_model.Domain.FromValues(P), f"b{x}") for x in cells}
for x in cells: m.Add(B[x] == sum(A[y] for y in nb[x]))
m.AddAllDifferent(list(B.values()))
m.Minimize(sum(B.values()))
s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = TL; s.parameters.num_workers = 2
st = s.Solve(m)
print(f"N={N}: status {s.StatusName(st)}, sum(B) = {s.ObjectiveValue():.0f}, lower bound {s.BestObjectiveBound():.0f}, {s.WallTime():.0f}s")
if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    for r in range(N): print('  A:', ' '.join(f"{s.Value(A[(r,c)]):3d}" for c in range(N)), '   B:', ' '.join(f"{s.Value(B[(r,c)]):4d}" for c in range(N)))
