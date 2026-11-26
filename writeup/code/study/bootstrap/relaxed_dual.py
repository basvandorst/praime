# Dual solution of the k=1 convex program in the Fourier basis: rho = sum_v y_v v v^T (y = LP duals of the cuts).
# Prints g(log x)/g(0) at chosen x, and the grid points where the certificate is closest to zero.
import sys, numpy as np

from scipy.optimize import linprog
from conductor import FormG
X = float(sys.argv[1]); kmax = int(sys.argv[2]); rounds = int(sys.argv[3])
N = max(16, int(10 + 4*np.log(X)/np.log(2)))
F = FormG(X, N, [0]); two_a = 2*F.a; A = F.A
grid = [1.02]
while grid[-1]*1.004 < X - 1e-6: grid.append(grid[-1]*1.004)
grid = np.array(sorted(set(grid) | {float(m) for m in range(2, int(X) + 1) if m < X - 1e-6}))
S1 = np.array([F.S(np.log(x)) for x in grid])            # g_v(log x) = v' S v (unnormalised)
Cs = np.empty_like(S1)
for j, x in enumerate(grid):
    M = np.zeros_like(A); lx = np.log(x); k = 1
    while k*lx <= two_a and k <= kmax: M += lx*x**(-k/2)*F.S(k*lx); k += 1
    Cs[j] = M
w, V = np.linalg.eigh(A)
cuts = [V[:, i] for i in range(V.shape[1])]
a_v = [float(v @ A @ v) for v in cuts]; c_v = [np.einsum('i,jik,k->j', v, Cs, v) for v in cuts]
for r in range(rounds):
    cost = np.zeros(len(grid) + 1); cost[-1] = -1
    Aub = np.hstack([np.array(c_v), two_a*np.ones((len(cuts), 1))])
    res = linprog(cost, A_ub=Aub, b_ub=np.array(a_v), bounds=[(0, 6.0)]*len(grid) + [(None, None)], method='highs')
    n = res.x[:-1]; t = res.x[-1]
    H = A - np.einsum('j,jik->ik', n, Cs); w, V = np.linalg.eigh(H); lam = w[0]/two_a
    if lam >= t - 2e-9: break
    for i in range(3):
        v = V[:, i]; cuts.append(v); a_v.append(float(v @ A @ v)); c_v.append(np.einsum('i,jik,k->j', v, Cs, v))
y = -res.ineqlin.marginals; y = np.maximum(y, 0)
rho = sum(yi*np.outer(v, v) for yi, v in zip(y, cuts) if yi > 0)
print('X', X, 'kmax', kmax, 'rounds', r, 'lam', lam, 'LP t', t, 'sum y*2a', y.sum()*two_a, 'rank(rho)', int((np.linalg.eigvalsh(rho) > 1e-9*np.trace(rho)).sum()))
print('lowest eigenvalues of H/2a:', w[:5]/two_a)
g0 = np.trace(rho)*two_a/ two_a   # S(0) = 2a * identity?  check below
S0 = F.S(0.0); g0 = float(np.sum(rho*S0))
gx = np.einsum('ik,jik->j', rho, S1)/g0
cx = np.einsum('ik,jik->j', rho, Cs)/g0
print('W(g)/g(0) =', float(np.sum(rho*A))/g0, '   min certificate on grid:', cx.min())
for x0 in [1.02, 1.1, 1.2, 1.33, 1.5, 1.75, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5]:
    if x0 < X:
        j = int(np.argmin(abs(grid - x0))); print('  x=%.3f  g(log x)/g(0) = % .3e   certificate = % .3e   weight %.3f' % (grid[j], gx[j], cx[j], n[j]))
loc = [j for j in range(1, len(grid)-1) if cx[j] <= cx[j-1] and cx[j] <= cx[j+1]]
print('local minima of the certificate (x, value):', ' '.join('%.3f:%.1e' % (grid[j], cx[j]) for j in loc))
