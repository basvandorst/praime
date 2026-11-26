# Theta family with B(t) = sum_k c_k cos(k delta t): g = g0 * beta, beta = sum c_k (delta_{a_k}+delta_{-a_k})/2.
import numpy as np, sys
from scipy.optimize import linprog
sys.path.insert(0, '/tmp/claude-1000/conv_prover'); from common import *
delta = float(sys.argv[1]); K = int(sys.argv[2])
a = delta*np.arange(K+1)
# objective o_k = cosh(a_k/2)/2 + (1/pi) int_0^inf Xi^2 cos(a_k t) Omega dt
o = np.cosh(a/2)/2 + (np.cos(np.outer(a, t)) @ (w*X2*Om))/np.pi
gnorm = g0(a)                       # g(0) = sum c_k g0(a_k)
U = a[-1] + 5.0
uu = np.arange(0, U, 0.004)
G0 = {}
def g0s(x): return g0(x)
# g(u) = sum_k c_k (g0(u-a_k)+g0(u+a_k))/2 ; precompute g0 on a fine grid and interpolate
grid = np.arange(-U-a[-1]-1, U+a[-1]+1, 0.002); g0grid = g0(np.abs(grid))
def g0i(x): return np.interp(x, grid, g0grid)
Mg = np.stack([(g0i(uu-ak)+g0i(uu+ak))/2 for ak in a], axis=1)
th = np.linspace(0, np.pi, 4000)
MB = np.cos(np.outer(th, np.arange(K+1)))
res = linprog(o, A_ub=np.vstack([-Mg, -MB]), b_ub=np.zeros(uu.size+th.size), A_eq=gnorm[None,:], b_eq=[1.0],
              bounds=[(None,None)]*(K+1), method='highs')
c = res.x
print("delta=%g K=%d A=%g status=%s  W/g(0) = %.6e" % (delta, K, a[-1], res.status, res.fun))
# check with prime sum
pp = lam_table(400)
gl = np.array([sum(c[k]*(g0i(np.log(n)-a[k])+g0i(np.log(n)+a[k]))/2 for k in range(K+1)) for n,_ in pp])
print(" prime side 2 sum Lambda n^-1/2 g(log n) =", 2*sum(L*n**-0.5*v for (n,L),v in zip(pp,gl)))
print(" min g on grid:", (Mg@c).min(), " min B:", (MB@c).min())
big = [(n, v) for (n,_),v in zip(pp,gl) if abs(v) > 1e-6]
print(" g(log n) > 1e-6 at:", [(n, round(v,6)) for n,v in big[:15]])
np.save('/tmp/claude-1000/conv_prover/lpcos_c_%g_%d.npy' % (delta,K), c)
