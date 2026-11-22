"""Direct time average over u = ln X of the off-diagonal part of R(X,d) (large-X form):
off(u) = (2/d) Re sum_{j!=k, both signs} c_j conj(c_k) e^{i(g_j-g_k)u}/(2+i(g_j-g_k)).
Compare mean square with (4/d^2)*Q and (4/d^2)*2Q. Also check against the finite-X r93 form at a few X."""
import numpy as np, math, sys
from q import Qparts
Z = np.loadtxt('/home/bas/praime/writeup/data/riemann/zeros_10000_flint.txt')
M = int(sys.argv[1]); NU = int(sys.argv[2])
g = Z[:M]; G = np.concatenate([g, -g]); rho = 0.5 + 1j*G
rng = np.random.default_rng(12345)
for d in (0.1, 0.03, 0.01):
    c = ((1+d)**rho - 1)/rho
    dg = G[:, None] - G[None, :]
    A = 1/(2 + 1j*dg); np.fill_diagonal(A, 0)
    u = rng.uniform(10, 10 + 1e6, NU)
    vals = []
    for i0 in range(0, NU, 1000):
        V = c[None, :]*np.exp(1j*np.outer(u[i0:i0+1000], G))
        vals.append(np.real(np.sum((V @ A)*np.conj(V), 1)))
    off = 2/d*np.concatenate(vals)
    _, same, opp, selfc = Qparts(g, ds=(d,))
    Q = 2*same[0] + 2*opp[0]
    ms = np.mean(off**2); se = np.std(off**2)/math.sqrt(NU)
    print('d', d, 'M', M, 'mean', np.mean(off), 'MS', ms, '+-', se, '| (4/d2)Q', 4/d**2*Q, '(4/d2)2Q', 8/d**2*Q,
          'LI exact', 4/d**2*(2*Q - 2*selfc[0]), '| rms', math.sqrt(ms), 'pred', 2/d*math.sqrt(2*Q), flush=True)
    # finite-X r93 form check of the off-diagonal at one u (X = e^u): include -1 term
    u0 = 30.0
    num = np.exp((2 + 1j*dg)*u0) - 1.0
    full = (np.conj(c)[None, :]*c[:, None]*num/(2+1j*dg)); np.fill_diagonal(full, 0)
    r93off = 2*np.real(full.sum())/(d*math.exp(2*u0))
    V = c*np.exp(1j*u0*G); mine = 2/d*np.real(V @ A @ np.conj(V))
    print('   r93-form offdiag at u=30', r93off, 'mine', mine)
