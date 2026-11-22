import numpy as np, math
from scipy.linalg import eigh_tridiagonal
Z = np.loadtxt('/home/bas/praime/writeup/data/riemann/zeros_10000_flint.txt')
def Nbar(t): return t/(2*math.pi)*np.log(t/(2*math.pi*math.e)) + 7/8
def Ninv(y):
    y = np.asarray(y, float)
    t = np.maximum(2*math.pi*np.exp(1.0) * np.ones_like(y), 2*math.pi*np.maximum(y, 1)/np.log(np.maximum(y, 3)) * 1.0 + 10)
    for _ in range(60):
        t = t - (Nbar(t) - y)/(np.log(t/(2*math.pi))/(2*math.pi))
        t = np.maximum(t, 2*math.pi*1.0001)
    return t
N1 = float(Nbar(Z[0])); NT = float(Nbar(Z[-1]))
def gue_unfolded(seed, n=30000, radius='exact'):
    rng = np.random.default_rng(seed)
    dgl = rng.normal(0, math.sqrt(2), n)
    off = np.sqrt(rng.chisquare(2*np.arange(n-1, 0, -1)))
    ev = np.sort(eigh_tridiagonal(dgl, off, eigvals_only=True))/math.sqrt(2*n)
    R = 2.0 if radius == 'exact' else float(np.max(np.abs(ev)))*(1+1e-9)
    x = np.clip(ev/R, -1, 1)
    return n*(0.5 + (x*np.sqrt(1-x*x) + np.arcsin(x))/math.pi), ev
