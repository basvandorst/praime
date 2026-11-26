# Convex relaxation: real weights n_j >= 0 on a grid of positions instead of a multiset of points.
# The form is linear in the weights, so max lambda_min is a concave problem: Kelley cutting planes with
# an LP (HiGHS).  The optimum is global (on the grid).  Float64.
import sys, numpy as np
from scipy.optimize import linprog
from conductor import FormG
def solve(X, kappas, step=0.01, cap=4.0, rounds=400, lo=1.02):
    N = max(12, int(8 + 4*np.log(X)/np.log(2)))
    F = FormG(X, N, kappas); two_a = 2*F.a
    grid = np.unique(np.round(np.concatenate([np.arange(lo, X - 1e-9, step), np.arange(2, int(X) + 1)]), 6))
    grid = grid[grid < X - 1e-9]
    C = np.array([F.comb_one(x) if hasattr(F, 'comb_one') else None for x in grid]) if False else None
    Cs = []
    for x in grid:
        M = np.zeros_like(F.A); lx = np.log(x); k = 1
        while k*lx <= 2*F.a:
            M += lx*x**(-k/2)*F.S(k*lx); k += 1
        Cs.append(M)
    Cs = np.array(Cs); A = F.A
    w, V = np.linalg.eigh(A); cuts = [V[:, i] for i in range(V.shape[1])]
    n = np.zeros(len(grid)); t = None
    for r in range(rounds):
        Vc = np.array(cuts)
        a_v = np.einsum('ci,ij,cj->c', Vc, A, Vc)
        c_vj = np.einsum('ci,jik,ck->cj', Vc, Cs, Vc)
        # maximise t  s.t.  a_v - sum_j n_j c_vj >= t * 2a
        cost = np.zeros(len(grid) + 1); cost[-1] = -1
        Aub = np.hstack([c_vj, two_a*np.ones((len(cuts), 1))]); bub = a_v
        res = linprog(cost, A_ub=Aub, b_ub=bub, bounds=[(0, cap)]*len(grid) + [(None, None)], method='highs')
        n = res.x[:-1]; t = res.x[-1]
        H = A - np.einsum('j,jik->ik', n, Cs)
        w, V = np.linalg.eigh(H); lam = w[0]/two_a
        if lam >= t - 1e-13: break
        for i in range(3): cuts.append(V[:, i])
    return grid, n, lam, t, r
if __name__ == '__main__':
    kappas = [int(c) for c in sys.argv[1]]
    for X in [float(v) for v in sys.argv[2].split(',')]:
        grid, n, lam, t, r = solve(X, kappas)
        nz = [(round(float(x), 3), round(float(v), 6)) for x, v in zip(grid, n) if v > 1e-7]
        print('kappas %s X=%g: max lambda_min over real weights = %.3e (LP bound %.3e, %d rounds); conductor %.8f' % (sys.argv[1], X, lam, t, r, np.exp(-lam)))
        print('    weights (position, weight):', nz, flush=True)
