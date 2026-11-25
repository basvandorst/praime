import numpy as np, sys
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
from analysis import *
def joint(a, rows=(1, 2, 3)):
    X = a['X']; B, A, Cres, _ = ols(X, a['k3']); pf = ols(X, a['full'])[0]; pq = ols(X, a['Q0'])[0]; pdg = ols(X, a['diagonly'])[0]
    idx = [(r, j) for r in rows for j in range(16)]
    V = np.array([[A[r, r2]*Cres[j, j2] for (r2, j2) in idx] for (r, j) in idx]); Vi = np.linalg.inv(V)
    dQ = np.array([pf[r, j] - pq[r, j] for (r, j) in idx]); dO = np.array([pf[r, j] - pdg[r, j] for (r, j) in idx])
    y = np.array([B[r, j] - pf[r, j] for (r, j) in idx]) + dQ + dO
    M = np.vstack([dQ, dO]).T; Cov = np.linalg.inv(M.T @ Vi @ M); g = Cov @ M.T @ Vi @ y
    e = y - M @ g
    return g, np.sqrt(np.diag(Cov)), Cov[0, 1]/np.sqrt(Cov[0, 0]*Cov[1, 1]), e @ Vi @ e, len(idx) - 2
Dsel = [D for D in Ds if D != 1153]; a = arrays(Dsel)
for rows, lab in (((1, 2, 3), 'p=2,3,5'), ((1,), 'p=2'), ((2,), 'p=3'), ((3,), 'p=5'), ((1, 2, 3, 4, 5, 6), 'p=2..13')):
    g, s, rho, chi, dof = joint(a, rows)
    print('%-8s joint GLS: gamma_Q %.3f +- %.3f, gamma_off %.3f +- %.3f, corr %.2f, chi2 %.1f / %d' % (lab, g[0], s[0], g[1], s[1], rho, chi, dof))
