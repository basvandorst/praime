import numpy as np, mpmath as mp, os
D = '/tmp/claude-1000/conv_prover/'
d = np.load(D+'xi_grid.npz'); t = d['t']; X = d['xi']; h = t[1]-t[0]
X2 = X**2
w = np.full(t.size, h); w[0] = h/2
if os.path.exists(D+'omega.npy'):
    Om = np.load(D+'omega.npy')
else:
    Om = np.array([float(mp.re(mp.digamma(mp.mpf(1)/4 + 1j*mp.mpf(x)/2))) for x in t]) - np.log(np.pi)
    np.save(D+'omega.npy', Om)
def ft_even(weights_t, u):
    """(1/2pi) int_R F(t) cos(tu) dt for even F given on the half grid (F*w precomputed)."""
    u = np.atleast_1d(np.asarray(u, float))
    out = np.empty(u.size)
    for i in range(0, u.size, 400):
        out[i:i+400] = (np.cos(np.outer(u[i:i+400], t)) @ weights_t)/np.pi
    return out
def g0(u):  return ft_even(w*X2, u)
def gk(u, k):  return ft_even(w*X2*t**(2*k), u)
import sympy
def lam_table(N):
    out = []
    for n in range(2, N+1):
        f = sympy.factorint(n)
        if len(f) == 1: out.append((n, np.log(list(f)[0])))
    return out
