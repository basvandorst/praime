# Puzzle 764 (writeup/764.md): a lower bound for large N by relaxing the path structure. Every label g (1..N-1) needs
# two distinct cells that are multiples of g (the two ends of its edge), each cell carries at most two labels, and the
# N cells are distinct positive integers <= B = S - N(N-1)/2. Integer program (HiGHS):
#   x_v = 1 if value v is a cell;  y_{g,v} = 1 if a copy of label g sits in cell v (only if g | v);
#   sum_v y_{g,v} = 2;  sum_g y_{g,v} <= 2 x_v;  sum_v x_v = N;  minimise sum v x_v.
# Usage: python3 p764_cover_bound.py N S [time_s]
import sys, numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix

N, S = int(sys.argv[1]), int(sys.argv[2]); tl = float(sys.argv[3]) if len(sys.argv) > 3 else 600
B = S - N * (N - 1) // 2; K = N - 1
pairs = [(g, v) for g in range(1, K + 1) for v in range(g, B + 1, g)]
nx = B; ny = len(pairs); n = nx + ny
print('N=%d S=%d: values <= %d, %d (label, value) pairs' % (N, S, B, ny), flush=True)
A = lil_matrix((K + B + 1, n)); lo, hi = [], []
r = 0
for g in range(1, K + 1):                       # each label: exactly two cells
    for j, (gg, v) in enumerate(pairs):
        if gg == g: A[r, nx + j] = 1
    lo.append(2); hi.append(2); r += 1
for v in range(1, B + 1):                       # each cell: at most two labels, only if used
    A[r, v - 1] = -2
    for j, (gg, vv) in enumerate(pairs):
        if vv == v: A[r, nx + j] = 1
    lo.append(-np.inf); hi.append(0); r += 1
for v in range(1, B + 1): A[r, v - 1] = 1        # N cells
lo.append(N); hi.append(N); r += 1
c = np.zeros(n); c[:nx] = np.arange(1, B + 1)
res = milp(c, constraints=LinearConstraint(A.tocsr(), lo, hi), integrality=np.ones(n), bounds=Bounds(0, 1),
           options={'time_limit': tl, 'mip_rel_gap': 1e-6})
print('status:', res.message)
print('cover bound: objective %.1f, solver lower bound %.1f  (record %d)' % (res.fun, res.mip_dual_bound, S))
