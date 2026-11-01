# Problem 68 (writeup/prob-68.md): second, independent solve of the minimum-P program of prob68_ilp.py with a
# different solver (OR-Tools CP-SAT) and separately written model code. Same reduction: k = 2 (mod 3), the +1 side
# must cover odd n = 2t + 1, the -1 side even n = 2t, t mod T; a prime of order m covers
#   m odd: one class mod m on one side;  4 | m: one class mod m/2 on one side;
#   m = 2 (mod 4), h = m/2: class c mod h on the +1 side and c + (h + 1)/2 mod h on the -1 side.
# A paired prime counts as shared when it has a private t on both sides. Cost: round(10^6 * log10 p) per prime.
# Usage: python3 prob68_cpsat.py T s [time_limit_s]    (T = 360, s = 1: about a minute; s = 2: a few minutes)
import sys, os, json, math
from ortools.sat.python import cp_model
from sympy import n_order

T, s_min = int(sys.argv[1]), int(sys.argv[2]); tl = float(sys.argv[3]) if len(sys.argv) > 3 else 3600
here = os.path.dirname(os.path.abspath(__file__))
orders = {int(m): ps for m, ps in json.load(open(os.path.join(here, '..', 'data', 'prob68_primes.json'))).items()}
M = cp_model.CpModel()
one, paired = {}, {}                  # one[p] = (d, {(side, c): var}); paired[p] = (h, {c: var})
for m, ps in orders.items():
    for p in ps:
        if p == 3: continue
        d = m if m % 2 else m // 2
        if T % d: continue
        if m % 4 == 2: paired[p] = (d, {c: M.NewBoolVar('') for c in range(d)})
        else: one[p] = (d, {(sd, c): M.NewBoolVar('') for sd in 'AB' for c in range(d)})
used = {}
for p, (d, vs) in list(one.items()) + list(paired.items()):
    u = M.NewBoolVar(''); M.Add(sum(vs.values()) == u); used[p] = u
def coverers(side, t):
    out = []
    for p, (d, vs) in one.items(): out.append(vs[(side, t % d)])
    for p, (h, vs) in paired.items():
        c = t % h if side == 'A' else (t - (h + 1) // 2) % h
        out.append(vs[c])
    return out
cov = {(sd, t): coverers(sd, t) for sd in 'AB' for t in range(T)}
for key, lits in cov.items(): M.AddBoolOr(lits)
shared = []
for p, (h, vs) in paired.items():
    sh = M.NewBoolVar(''); shared.append(sh)
    for sd in 'AB':
        priv = []
        for t in range(T):
            mine = vs[t % h] if sd == 'A' else vs[(t - (h + 1) // 2) % h]
            z = M.NewBoolVar(''); priv.append(z)
            M.AddImplication(z, mine)
            others = [l for l in cov[(sd, t)] if l is not mine]
            for l in others: M.AddImplication(z, l.Not())
        M.AddBoolOr(priv).OnlyEnforceIf(sh)
M.Add(sum(shared) >= s_min - 1)
for p in (5, 7):                       # translations t -> t + 1, t -> t + 2: fix the classes of 5 and 7 to 0
    if p in one:
        for (sd, c), v in one[p][1].items():
            if c: M.Add(v == 0)
M.Minimize(sum(round(1e6 * math.log10(p)) * u for p, u in used.items()))
solver = cp_model.CpSolver(); solver.parameters.max_time_in_seconds = tl; solver.parameters.num_workers = 2
st = solver.Solve(M)
chosen = sorted([3] + [p for p, u in used.items() if solver.Value(u)])
P = math.prod(chosen)
print('T = %d, s >= %d: %s; P = %d (%d digits), primes %s' % (T, s_min, solver.StatusName(st), P, len(str(P)), chosen))
print('objective %.6f, best bound %.6f (log10 P minus log10 3, in units of 1e-6)' % (solver.ObjectiveValue() / 1e6, solver.BestObjectiveBound() / 1e6))
