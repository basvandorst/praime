# Puzzle 772, Q2: most primes in A and B together (parity bound), then minimal sum(B), one CP-SAT run per parity
# configuration that reaches the bound. Usage: p772_q2_cpsat.py N U TIMELIMIT   (N = 3: OPTIMAL 514)
import sys
from ortools.sat.python import cp_model
from sympy import primerange
N, U, TL = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
cells = [(r, c) for r in range(N) for c in range(N)]; n = len(cells)
nb = [[j for j, y in enumerate(cells) if abs(x[0]-y[0]) <= 1 and abs(x[1]-y[1]) <= 1] for x in cells]
w = [sum(1 for c in range(n) if j in nb[c]) for j in range(n)]
def total(E):
    ne = bin(E).count('1'); pa = (n - ne) + (1 if ne else 0)
    ob = [c for c in range(n) if sum(1 for j in nb[c] if not (E >> j & 1)) % 2]
    return pa + len(ob), ob
mx = max(total(E)[0] for E in range(1 << n))
configs = [E for E in range(1 << n) if total(E)[0] == mx]
print(f"N={N}: max primes {mx}, {len(configs)} parity configurations")
oddP = list(primerange(3, U + 1)); evens = list(range(0, U + 1, 2)); P = set(primerange(2, 9 * U + 10))
best = None
for E in configs:
    evc = [j for j in range(n) if E >> j & 1]; ob = total(E)[1]
    for two in (evc if evc else [None]):
        m = cp_model.CpModel()
        A = []
        for j in range(n):
            if E >> j & 1:
                A.append(m.NewConstant(2) if j == two else m.NewIntVarFromDomain(cp_model.Domain.FromValues([e for e in evens if e != 2]), f"a{j}"))
            else:
                A.append(m.NewIntVarFromDomain(cp_model.Domain.FromValues(oddP), f"a{j}"))
        m.AddAllDifferent(A)
        B = []
        for c in range(n):
            if c in ob: b = m.NewIntVarFromDomain(cp_model.Domain.FromValues(sorted(P - {2})), f"b{c}")
            else: b = m.NewIntVar(0, 9 * U, f"b{c}")
            m.Add(b == sum(A[j] for j in nb[c])); B.append(b)
        m.AddAllDifferent(B)
        m.Minimize(sum(B))
        s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = TL; s.parameters.num_workers = 2
        st = s.Solve(m)
        if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            v = s.ObjectiveValue()
            print(f"  even cells {[cells[j] for j in evc]}, 2 at {cells[two] if two is not None else None}: sum {v:.0f} ({s.StatusName(st)}, bound {s.BestObjectiveBound():.0f})")
            if best is None or v < best[0]: best = (v, [s.Value(a) for a in A], [s.Value(b) for b in B], s.StatusName(st))
        else:
            print(f"  even cells {[cells[j] for j in evc]}, 2 at {cells[two] if two is not None else None}: {s.StatusName(st)}")
v, A, B, st = best
print(f"best sum {v:.0f}:"); 
for r in range(N): print('  A:', ' '.join(f"{A[r*N+c]:3d}" for c in range(N)), '   B:', ' '.join(f"{B[r*N+c]:4d}" for c in range(N)))
