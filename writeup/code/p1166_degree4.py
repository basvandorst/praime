"""Puzzle 1166 (Nomoto): sets X of even numbers with  prod_{x in X} P_n(x) = (p_n# / 2)^4.

P_n(x) is the product of p (p + x) over all pairs of odd primes p < p + x <= p_n.  A prime q
appears in P_n(x) once for each of q - x and q + x that is an odd prime <= p_n, so the product
over X equals (p_n#/2)^4 exactly when every odd prime q <= p_n has exactly four 'partners'
q +- x (x in X) among the odd primes <= p_n.  That is an exact 0/1 feasibility problem:
   for every odd prime q <= p_n:  sum_{x in X} ([q - x odd prime] + [q + x prime <= p_n]) = 4.
We solve it with OR-Tools CP-SAT for each n and enumerate all solutions.
Usage: p1166_degree4.py NMAX [DEG]
"""
import sys
from sympy import prime, primerange
from ortools.sat.python import cp_model

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 40
DEG = int(sys.argv[2]) if len(sys.argv) > 2 else 4

NMIN = int(sys.argv[3]) if len(sys.argv) > 3 else 3
for n in range(NMIN, NMAX + 1):
    pn = prime(n)
    odd = [q for q in primerange(3, pn + 1)]
    S = set(odd)
    xs = list(range(2, pn - 3 + 1, 2))
    m = cp_model.CpModel()
    X = {x: m.NewBoolVar(f"x{x}") for x in xs}
    for q in odd:
        terms = []
        for x in xs:
            c = (q - x in S) + (q + x in S)
            if c:
                terms.append(c * X[x])
        m.Add(sum(terms) == DEG)
    solver = cp_model.CpSolver()
    solver.parameters.num_workers = 1
    solver.parameters.enumerate_all_solutions = True
    solver.parameters.max_time_in_seconds = 120

    class Collect(cp_model.CpSolverSolutionCallback):
        def __init__(self):
            super().__init__(); self.sols = []
        def on_solution_callback(self):
            self.sols.append([x for x in xs if self.Value(X[x])])
    cb = Collect()
    st = solver.Solve(m, cb)
    name = solver.StatusName(st)
    print(f"n = {n:3d}  p_n = {pn:4d}  status {name:10s}  solutions {len(cb.sols)}" +
          ("  " + "; ".join(",".join(map(str, s)) for s in cb.sols[:6]) if cb.sols else ""), flush=True)
