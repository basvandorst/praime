# Problem 68 (writeup/prob-68.md): smallest P = prod(A u B) over Brier numbers whose minimal covering sets A (for
# k*2^n + 1) and B (for k*2^n - 1) share at least s primes, 3 among them, as an exact integer program (HiGHS via scipy).
#
# Model. Take k = 2 (mod 3): then 3 divides k*2^n + 1 for even n and k*2^n - 1 for odd n. What is left: A must cover
# the odd n = 2t + 1 and B the even n = 2t, for t mod T (so n mod 2T). For a prime p with m = ord_p(2):
#   m odd        : one class mod m, on the A side or on the B side (k in -<2> or in <2> mod p);
#   m = 0 mod 4  : one class mod m/2 on one side (its class on the other side has the parity that 3 already covers);
#   m = 2 mod 4  : "paired": class c mod h = m/2 on the A side and c + 1/2 (mod h) on the B side, always both.
# Only paired primes (and 3) can be in A and B at once. A prime is in a minimal covering set iff it covers some point
# that no other member covers (a private point); the program asks for >= s - 1 paired primes with private points on
# both sides and minimises sum log10 p.
#
# Usage: python3 prob68_ilp.py T s [time_limit_s] [shared_pool] [excl]
#   T           t-modulus (360 for periods n | 720; 1260 for n | 2520)
#   s           required |A n B|, 3 included
#   shared_pool comma list of paired primes that may count as shared (default: all paired primes)
#   excl        drop the paired primes outside the pool altogether (then |A n B| <= 1 + |pool|)
# Prints the optimum (or best found within the time limit), the solver's lower bound, and k, A, B built by CRT.
import sys, os, json, math, numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix
from sympy import n_order
from sympy.ntheory.modular import crt

T = int(sys.argv[1]); s_min = int(sys.argv[2])
tl = float(sys.argv[3]) if len(sys.argv) > 3 else 1500
pool = [int(x) for x in sys.argv[4].split(',')] if len(sys.argv) > 4 else None
excl = 'excl' in sys.argv
here = os.path.dirname(os.path.abspath(__file__))
pr = json.load(open(os.path.join(here, '..', 'data', 'prob68_primes.json')))

cands = []                                   # (p, kind, d)
for m, ps in pr.items():
    m = int(m)
    for p in ps:
        if p == 3: continue
        assert n_order(2, p) == m
        kind, d = ('one', m) if m % 2 else (('one', m // 2) if m % 4 == 0 else ('pair', m // 2))
        if T % d: continue
        if kind == 'pair' and excl and pool is not None and p not in pool: continue
        cands.append((p, kind, d))
cands.sort()
var = []                                     # (p, kind, d, side, c), side 'A', 'B' or 'P' (paired)
for p, kind, d in cands:
    for c in range(d):
        var += [(p, kind, d, 'A', c), (p, kind, d, 'B', c)] if kind == 'one' else [(p, kind, d, 'P', c)]
nx = len(var)
pairs = [p for p, kind, d in cands if kind == 'pair' and (pool is None or p in pool)]

def covers(v, side, t):
    p, kind, d, sd, c = v
    if sd == 'P': return t % d == (c if side == 'A' else (c + (d + 1) // 2) % d)
    return sd == side and t % d == c

zidx = {}
for p in pairs:
    for side in 'AB':
        for t in range(T): zidx[(p, side, t)] = nx + len(zidx)
shidx = {p: nx + len(zidx) + i for i, p in enumerate(pairs)}
N = nx + len(zidx) + len(shidx)
cov = {(side, t): [] for side in 'AB' for t in range(T)}
byp = {}
for i, v in enumerate(var):
    byp.setdefault(v[0], []).append(i)
    for side in 'AB':
        if v[3] in ('P', side):
            for t in range(v[4] % v[2], T, v[2]) if v[3] != 'P' or side == 'A' else \
                     range((v[4] + (v[2] + 1) // 2) % v[2], T, v[2]):
                cov[(side, t)].append(i)
M = lil_matrix((len(byp) + 2 * T + 2 * len(zidx) + len(pairs) + 20, N)); lo, hi = [], []; r = 0
for p, idx in byp.items():                   # one class per prime
    for i in idx: M[r, i] = 1
    lo.append(-np.inf); hi.append(1); r += 1
for key, idx in cov.items():                 # both sides covered
    for i in idx: M[r, i] = 1
    lo.append(1); hi.append(np.inf); r += 1
BIG = len(cands) + 1
for (p, side, t), zi in zidx.items():        # z = 1: t is a private point of p on this side
    for i in byp[p]:
        if covers(var[i], side, t): M[r, i] -= 1
    M[r, zi] = 1; lo.append(-np.inf); hi.append(0); r += 1
    for i in cov[(side, t)]: M[r, i] = 1
    M[r, zi] = BIG; lo.append(-np.inf); hi.append(BIG + 1); r += 1
for p in pairs:                              # shared: a private point on each side
    for side in 'AB':
        for t in range(T): M[r, zidx[(p, side, t)]] = 1
        M[r, shidx[p]] = -1; lo.append(0); hi.append(np.inf); r += 1
for p in pairs: M[r, shidx[p]] = 1
lo.append(s_min - 1); hi.append(np.inf); r += 1
for i, v in enumerate(var):                  # translations t -> t + j: fix the classes of 5 and 7 to 0
    if v[0] in (5, 7) and v[4] != 0:
        M[r, i] = 1; lo.append(0); hi.append(0); r += 1
cost = np.zeros(N)
for i, v in enumerate(var): cost[i] = math.log10(v[0])
res = milp(cost, constraints=LinearConstraint(M[:r].tocsr(), lo, hi), integrality=np.ones(N), bounds=Bounds(0, 1),
           options={'time_limit': tl, 'mip_rel_gap': 0})
print('T = %d, s >= %d: %s; lower bound %.3f digits' % (T, s_min, res.message, math.log10(3) + res.mip_dual_bound))
if res.x is None: sys.exit()
sol = [var[i] for i in range(nx) if res.x[i] > 0.5]
mods, rs = [3], [2]
for p, kind, d, side, c in sol:              # back to k: A side n = 2c+1 -> k = -2^-n; B side n = 2c -> k = 2^-n
    m = n_order(2, p)
    rs.append((-pow(2, -(2 * c + 1) % m, p)) % p if side in 'AP' else pow(2, -(2 * c) % m, p)); mods.append(p)
k, P = (int(x) for x in crt(mods, rs))
k = k if k % 2 else k + P
A = sorted([3] + [v[0] for v in sol if v[3] in 'AP']); B = sorted([3] + [v[0] for v in sol if v[3] in 'BP'])
print('P = %d (%.3f digits), %d primes' % (P, math.log10(P), len(mods)))
print('k =', k)
print('A candidates =', A); print('B candidates =', B)
print('(paired primes without a private point on one side drop out of that set; prob68_verify.py decides)')
