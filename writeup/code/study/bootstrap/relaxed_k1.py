# Variant of relaxed_fields.py: positive measure d psi with NO power structure (k = 1 term only).
# Dual certificates are then g >= 0 with ghat >= 0.  Question: is this larger class also sharp at conductor 1?
import sys, time, numpy as np

from scipy.optimize import linprog
from conductor import FormG
def solve(X, kappas, kmax, ratio=1.004, cap=6.0, rounds=400, tol=2e-8):
    N = max(16, int(10 + 4*np.log(X)/np.log(2)))
    F = FormG(X, N, kappas); two_a = 2*F.a; A = F.A
    grid = [1.02]
    while grid[-1]*ratio < X - 1e-6: grid.append(grid[-1]*ratio)
    grid = np.array(sorted(set(grid) | {float(m) for m in range(2, int(X) + 1) if m < X - 1e-6}))
    Cs = np.empty((len(grid),) + A.shape)
    for j, x in enumerate(grid):
        M = np.zeros_like(A); lx = np.log(x); k = 1
        while k*lx <= two_a and k <= kmax:
            M += lx*x**(-k/2)*F.S(k*lx); k += 1
        Cs[j] = M
    w, V = np.linalg.eigh(A)
    cuts = [V[:, i] for i in range(V.shape[1])]
    a_v = [float(v @ A @ v) for v in cuts]; c_v = [np.einsum('i,jik,k->j', v, Cs, v) for v in cuts]
    for r in range(rounds):
        cost = np.zeros(len(grid) + 1); cost[-1] = -1
        Aub = np.hstack([np.array(c_v), two_a*np.ones((len(cuts), 1))])
        res = linprog(cost, A_ub=Aub, b_ub=np.array(a_v), bounds=[(0, cap)]*len(grid) + [(None, None)], method='highs')
        n = res.x[:-1]; t = res.x[-1]
        H = A - np.einsum('j,jik->ik', n, Cs)
        w, V = np.linalg.eigh(H); lam = w[0]/two_a
        if lam >= t - tol: break
        for i in range(3):
            v = V[:, i]; cuts.append(v); a_v.append(float(v @ A @ v)); c_v.append(np.einsum('i,jik,k->j', v, Cs, v))
    return grid, n, lam, t, r
def clusters(grid, n):
    out = []; cur = None
    for x, v in zip(grid, n):
        if v > 1e-6:
            if cur is not None and x/cur[2] < 1.012: cur = [cur[0] + v*x, cur[1] + v, x]
            else:
                if cur is not None: out.append(cur)
                cur = [v*x, v, x]
    if cur is not None: out.append(cur)
    return ['%.3f(%.3f)' % (c[0]/c[1], c[1]) for c in out]
for X in [float(v) for v in sys.argv[2].split(',')]:
    t0 = time.time()
    grid, n, lam, t, r = solve(X, [0], int(sys.argv[1]))
    print('kmax %s X=%-4g lambda %.4e (LP bound %.4e, %d rounds, %.0fs) %s' % (sys.argv[1], X, lam, t, r, time.time()-t0, ' '.join(clusters(grid, n))), flush=True)
