"""Puzzle 772: checks the grids and the parity bound, independently of the CP-SAT scripts.

1. Our grids for N = 3 (Q1: all B prime, sum 291; Q2: 17 primes, sum 514) and the page's grids: A distinct
   non-negative, B = sum over the king neighbourhood including the cell, B distinct, prime counts and sums.
2. Parity bound for Q2: an even A value is prime only if it is 2, and B_c is odd iff an odd number of the A values
   around c are odd; maximising over all sets of even cells gives 17 (N = 3) and 29 (N = 4).
3. Optimality of the sums for N = 3 by an exhaustive argument independent of CP-SAT:
   Q1: sum(B) = 4*corners + 6*edges + 9*centre; every value <= 291/4, so a search over all A with values < 73.
"""
from itertools import permutations, combinations
from sympy import isprime

def B_of(A):
    N = len(A)
    return [[sum(A[r][c] for r in range(i-1, i+2) for c in range(j-1, j+2) if 0 <= r < N and 0 <= c < N) for j in range(N)] for i in range(N)]
def check(A, name):
    N = len(A); flat = [x for row in A for x in row]; B = B_of(A); fb = [x for row in B for x in row]
    assert min(flat) >= 0 and len(set(flat)) == N * N and len(set(fb)) == N * N
    pa = sum(isprime(x) for x in flat); pb = sum(isprime(x) for x in fb)
    print(f"{name}: sum(B) = {sum(fb)}, primes in A = {pa}, primes in B = {pb}")
    return sum(fb), pa, pb
assert check([[12, 0, 6], [14, 3, 2], [10, 4, 8]], "ours N=3 Q1") == (291, 2, 9)
assert check([[19, 3, 7], [5, 2, 11], [13, 17, 31]], "ours N=3 Q2") == (514, 9, 8)
assert check([[0, 4, 10], [6, 1, 2], [12, 18, 8]], "page N=3 Q1 (Vantieghem)")[0] == 309
assert check([[5, 9, 13], [3, 0, 1], [15, 11, 7]], "page N=3 (Rivera)") == (304, 5, 8)
assert check([[7, 6, 12, 5], [4, 0, 14, 22], [10, 2, 16, 20], [1, 24, 8, 3]], "page N=4 Q1 (Vantieghem)")[0] == 988
# 2. parity bound
for N, want in ((3, 17), (4, 29)):
    cells = [(r, c) for r in range(N) for c in range(N)]; n = len(cells)
    nb = [[j for j, y in enumerate(cells) if abs(x[0]-y[0]) <= 1 and abs(x[1]-y[1]) <= 1] for x in cells]
    best = 0
    for E in range(1 << n):
        ne = bin(E).count('1')
        best = max(best, (n - ne) + (1 if ne else 0) + sum(1 for c in range(n) if sum(1 for j in nb[c] if not (E >> j & 1)) % 2))
    assert best == want
    print(f"parity bound: at most {best} primes in A and B together for N = {N}")
# 3. exhaustive check of the N = 3 optima.  Q1: B all prime forces A's centre odd and the 8 others even (unique
#    parity solution); enumerate centre (odd), edges and corners (distinct evens) with 4*corners+6*edges+9*centre <= 291.
def best_q1(limit):
    best = None; ev = list(range(0, 74, 2))
    for ctr in range(1, 74, 2):
        if 9 * ctr > limit: break
        for E in combinations(ev, 4):
            se = sum(E)
            if 9 * ctr + 6 * se > limit: continue
            for C in combinations([e for e in ev if e not in E], 4):
                if 9 * ctr + 6 * se + 4 * sum(C) > limit: continue
                for pe in permutations(E):
                    for pc in permutations(C):
                        A = [[pc[0], pe[0], pc[1]], [pe[1], ctr, pe[2]], [pc[2], pe[3], pc[3]]]
                        B = [x for row in B_of(A) for x in row]
                        if len(set(B)) == 9 and all(isprime(x) for x in B):
                            s = sum(B)
                            if best is None or s < best: best = s
    return best
b = best_q1(291)
assert b == 291
print("exhaustive: no grid with all of B prime has sum below 291 (N = 3)")
print("all checks passed")

# 4. Q2 for N = 3 with a second solver (HiGHS via scipy.milp): centre = 2 (the only parity configuration reaching
#    17 primes), the other 8 cells distinct odd primes, the 8 non-centre B values distinct primes; minimise sum(B).
#    Each A cell picks one odd prime (binary per value), each non-centre B cell one prime; B = neighbourhood sum.
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from sympy import primerange
cells = [(r, c) for r in range(3) for c in range(3)]; ctr = (1, 1)
oc = [x for x in cells if x != ctr]
Ap = list(primerange(3, 130)); Bp = list(primerange(3, 1200))
nA, nB = len(oc) * len(Ap), len(oc) * len(Bp)
iA = lambda k, j: k * len(Ap) + j; iB = lambda k, j: nA + k * len(Bp) + j
nv = nA + nB; rows, lo, hi = [], [], []
def row(entries, l, h):
    r = np.zeros(nv)
    for i, v in entries: r[i] += v
    rows.append(r); lo.append(l); hi.append(h)
for k in range(len(oc)):
    row([(iA(k, j), 1) for j in range(len(Ap))], 1, 1)          # one value per A cell
    row([(iB(k, j), 1) for j in range(len(Bp))], 1, 1)          # one value per B cell
for j in range(len(Ap)): row([(iA(k, j), 1) for k in range(len(oc))], 0, 1)   # A distinct
for j in range(len(Bp)): row([(iB(k, j), 1) for k in range(len(oc))], 0, 1)   # B distinct
for k, x in enumerate(oc):                                       # B_x = 2 + sum of A over the neighbourhood
    ent = [(iB(k, j), Bp[j]) for j in range(len(Bp))]
    for kk, y in enumerate(oc):
        if abs(x[0]-y[0]) <= 1 and abs(x[1]-y[1]) <= 1: ent += [(iA(kk, j), -Ap[j]) for j in range(len(Ap))]
    row(ent, 2, 2)
cost = np.zeros(nv)
for k in range(len(oc)):
    for j in range(len(Bp)): cost[iB(k, j)] = Bp[j]
# objective: sum(B) = sum over non-centre B + centre B; centre B = 2 + sum of the 8 A values
for k in range(len(oc)):
    for j in range(len(Ap)): cost[iA(k, j)] += Ap[j]
res = milp(cost, constraints=[LinearConstraint(np.array(rows), lo, hi)], integrality=np.ones(nv), bounds=Bounds(0, 1),
           options={"time_limit": 600})
assert res.status == 0
print(f"HiGHS: minimal sum(B) for 17 primes, N = 3: {res.fun + 2:.0f}")
assert round(res.fun + 2) == 514
print("all checks passed (HiGHS)")
