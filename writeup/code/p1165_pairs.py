"""Puzzle 1165 (Nomoto): disjoint sets d, k of even distances with
    prod_{x in d} P_n(x) = prod_{x in k} P_n(x) = (p_n#/2)^m.
As in puzzle 1166, v_q(prod_{x in X} P_n(x)) = #{x in X: q - x odd prime} + #{x in X: q + x prime <= p_n}, so the
equation says: every odd prime q <= p_n has exactly m partners at distances in d, and exactly m at distances in k.
For each n and m this is a 0/1 feasibility problem; CP-SAT decides it (statuses printed), HiGHS rechecks n <= NH.
Usage: p1165_pairs.py NMAX NH TIMELIMIT
"""
import sys
import numpy as np
from sympy import prime, primerange
from ortools.sat.python import cp_model
from scipy.optimize import milp, LinearConstraint, Bounds

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 30
NH = int(sys.argv[2]) if len(sys.argv) > 2 else 22
TL = float(sys.argv[3]) if len(sys.argv) > 3 else 20

def system(n):
    pn = prime(n); odd = list(primerange(3, pn + 1)); S = set(odd)
    xs = [x for x in range(2, pn - 2, 2) if any(q + x in S for q in odd)]   # distances with at least one pair
    A = np.array([[(q - x in S) + (q + x in S) for x in xs] for q in odd])
    return pn, odd, xs, A

def cpsat(A, xs, m, enum=False):
    md = cp_model.CpModel(); X = len(xs)
    U = [md.NewBoolVar(f"u{i}") for i in range(X)]; V = [md.NewBoolVar(f"v{i}") for i in range(X)]
    for i in range(X): md.Add(U[i] + V[i] <= 1)
    for row in A:
        md.Add(sum(int(c) * U[i] for i, c in enumerate(row) if c) == m)
        md.Add(sum(int(c) * V[i] for i, c in enumerate(row) if c) == m)
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = TL
    if not enum:
        s.parameters.num_workers = 2; st = s.Solve(md)
        sol = ([xs[i] for i in range(X) if s.Value(U[i])], [xs[i] for i in range(X) if s.Value(V[i])]) if st in (cp_model.OPTIMAL, cp_model.FEASIBLE) else None
        return s.StatusName(st), sol
    s.parameters.num_workers = 1; s.parameters.enumerate_all_solutions = True
    class C(cp_model.CpSolverSolutionCallback):
        def __init__(s2): super().__init__(); s2.sols = []
        def on_solution_callback(s2): s2.sols.append((tuple(xs[i] for i in range(X) if s2.Value(U[i])), tuple(xs[i] for i in range(X) if s2.Value(V[i]))))
    cb = C(); st = s.Solve(md, cb); return s.StatusName(st), cb.sols

def highs(A, m):
    q, X = A.shape
    Aeq = np.zeros((2 * q, 2 * X)); Aeq[:q, :X] = A; Aeq[q:, X:] = A
    Adis = np.hstack([np.eye(X), np.eye(X)])
    cons = [LinearConstraint(Aeq, m, m), LinearConstraint(Adis, 0, 1)]
    r = milp(np.zeros(2 * X), constraints=cons, integrality=np.ones(2 * X), bounds=Bounds(0, 1), options={"time_limit": 120})
    return r.status            # 0 feasible (optimal), 2 infeasible

found = []
for n in range(3, NMAX + 1):
    pn, odd, xs, A = system(n)
    status = {}
    for m in range(1, len(odd)):
        st, sol = cpsat(A, xs, m)
        status[m] = st
        if sol: found.append((n, m, sol))
        if n <= NH:
            h = highs(A, m)
            assert (h == 0) == (sol is not None) and h in (0, 2), (n, m, st, h)
    bad = [m for m, st in status.items() if st not in ("INFEASIBLE", "OPTIMAL", "FEASIBLE")]
    print(f"n = {n} (p_n = {pn}): pairs for m = {[m for m in status if status[m] in ('OPTIMAL', 'FEASIBLE')]}; undecided m = {bad}", flush=True)
print("solutions:", [(n, m) for n, m, _ in found])
for n, m, _ in found:
    pn, odd, xs, A = system(n)
    st, sols = cpsat(A, xs, m, enum=True)
    print(f"n = {n}, m = {m}: {len(sols)} ordered pairs (d, k) (status {st}); e.g. d = {list(sols[0][0])}, k = {list(sols[0][1])}")
