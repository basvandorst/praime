# The dual certificate of the character problem (exploration).  The LP's dual weights y_c on the cuts v_c define a
# test-function mixture with g(w) = sum_c y_c g_{v_c}(w); printed: G(n) = sum_c y_c v_c^T S(log n) v_c = 2 g(log n) for the
# integers n <= X, next to the character.  If -chi(n) G(n) >= 0 at every prime power coprime to the conductor (and G = 0
# at the first power of the ramified prime), the local minimum over x in [-1, 1] sits at x = chi(p): the bound is attained.
# Usage: dual.py <kappa> <X>
import sys, numpy as np
from scipy.optimize import linprog
from sympy import kronecker_symbol
from maass import Problem
kappa = int(sys.argv[1]); X = float(sys.argv[2]); D = -3 if kappa else 5
grid = np.unique(np.concatenate([np.linspace(-1, 1, 81), np.linspace(-1, -0.98, 21), np.linspace(0.98, 1, 21), np.linspace(-0.02, 0.02, 21)]))
cf = lambda ph: ph; cf.pm = lambda m, x: x**m
P = Problem(0.0, 0, X, T=[('chi', [kappa], 1, cf)], th=grid, margin=14.0)
# rerun the cutting planes, keeping the cuts
nP, J = len(P.primes), len(grid); two_a = 2*P.B.a; rows, rhs, vs = [], [], []
w, V = np.linalg.eigh(P.A[0])
for i in range(12): a_v, row = P.cut(0, V[:, i]); rhs.append(a_v); rows.append(row); vs.append(V[:, i])
Aeq = np.zeros((nP, nP*J + 1))
for i in range(nP): Aeq[i, i*J:(i + 1)*J] = 1
cost = np.zeros(nP*J + 1); cost[-1] = -1
for rd in range(300):
    res = linprog(cost, A_ub=np.hstack([np.array(rows), two_a*np.ones((len(rows), 1))]), b_ub=np.array(rhs), A_eq=Aeq, b_eq=np.ones(nP), bounds=[(0, 1)]*(nP*J) + [(-50, 50)], method='highs')
    wts = res.x[:-1].reshape(nP, J); t = res.x[-1]; ev, V = np.linalg.eigh(P.forms(wts)[0]); added = 0
    for i in range(4):
        if ev[i]/two_a < t - 1e-8: a_v, row = P.cut(0, V[:, i]); rhs.append(a_v); rows.append(row); vs.append(V[:, i]); added += 1
    if not added: break
y = -res.ineqlin.marginals; y = y/y.sum()
print('kappa %d X %g: conductor %.8f, %d cuts, %d with dual weight > 1e-9' % (kappa, X, np.exp(-t), len(rows), int((y > 1e-9).sum())))
g0 = sum(yc*float(v @ P.B.S(0.0)[0] @ v) for yc, v in zip(y, vs))
print(' n  chi(n)  G(n)/G(1)   [prime powers marked *]')
from sympy import factorint
bad = 0
for n in range(2, int(X) + 1):
    G = sum(yc*float(v @ P.B.S(np.log(n))[0] @ v) for yc, v in zip(y, vs))/g0; ch = int(kronecker_symbol(D, n)); pp = len(factorint(n)) == 1
    flag = ''
    if pp and ch != 0 and -ch*G < -1e-7: flag = '  <-- sign condition fails'; bad += 1
    print('%3d  %+d  %+.6f %s%s' % (n, ch, G, '*' if pp else '', flag))
print('prime powers coprime to the conductor with -chi(n) G(n) < 0: %d' % bad)
