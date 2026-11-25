import numpy as np, sys
from scipy import stats
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
from analysis import *
Dsel = [D for D in Ds if D != 1153]; a = arrays(Dsel); n = len(Dsel); X = a['X']
_, A, Cy, _ = ols(X, a['k3']); _, _, Cf, _ = ols(X, a['full']); _, _, Cd, _ = ols(X, a['k3'] - a['full'])
print('residual sd x1e4: data k3     ', np.round(np.sqrt(np.diag(Cy))*1e4).astype(int).tolist())
print('                  law full    ', np.round(np.sqrt(np.diag(Cf))*1e4).astype(int).tolist())
print('                  k3 - full   ', np.round(np.sqrt(np.diag(Cd))*1e4).astype(int).tolist())
B = ols(X, a['k3'])[0]; pf = ols(X, a['full'])[0]; pq = ols(X, a['Q0'])[0]; pdg = ols(X, a['diagonly'])[0]
dof = n - 7
for r, p in zip(range(1, 7), PR):
    e = B[r] - pf[r]
    c1 = e @ np.linalg.solve(A[r, r]*Cy, e); c2 = e @ np.linalg.solve(A[r, r]*Cd, e)
    pv = lambda c: 1 - stats.f.cdf(c*(dof - 15)/(dof*16), 16, dof - 15)
    print('p=%2d chi2 vs law: registered cov %.1f (p %.3f) | cov of (k3-full) residuals %.1f (p %.3f)' % (p, c1, pv(c1), c2, pv(c2)))
# gamma_Q with the (k3 - full) residual covariance, GLS
idx = [(r, j) for r in (1, 2, 3) for j in range(16)]
for lab, Cc in (('registered Cres', Cy), ('Cres of k3-full', Cd)):
    V = np.array([[A[r, r2]*Cc[j, j2] for (r2, j2) in idx] for (r, j) in idx]); Vi = np.linalg.inv(V)
    for nm, palt in (('gQ', pq), ('goff', pdg)):
        dQ = np.array([pf[r, j] - palt[r, j] for (r, j) in idx]); obs = np.array([B[r, j] - palt[r, j] for (r, j) in idx])
        print('  %s GLS %s = %.4f +- %.4f' % (lab, nm, dQ @ Vi @ obs/(dQ @ Vi @ dQ), 1/np.sqrt(dQ @ Vi @ dQ)))
