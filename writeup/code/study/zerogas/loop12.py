import sys, numpy as np
from sympy import primerange
from loop9 import cth
from loop11 import cn, al
lo, hi, umax = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
g = np.loadtxt('/home/bas/praime/writeup/data/riemann/odlyzko_zeros1.txt')[lo:hi]
theta = lambda t: t/2*np.log(t/(2*np.pi)) - t/2 - np.pi/8 + 1/(48*t) + 7/(5760*t**3)
Sz = np.arange(lo + 1, hi + 1) - 0.5 - (theta(g)/np.pi + 1); Sz = Sz - Sz.mean(); L = float(np.log(g/(2*np.pi)).mean())
nmax = int(np.exp(umax*L)) + 1; terms = []
for p in primerange(2, nmax + 1):
    m = 1
    while p**m <= nmax: terms.append((p, m)); m += 1
terms.sort(key=lambda pm: pm[1]*np.log(pm[0]))
edges = np.arange(0, umax + 1e-9, 0.1); res = Sz.copy(); cum_th = 0.0; v0 = Sz.var()
print('zeros %d..%d: L = %.3f, Var S at the zeros = %.5f, %d prime powers to u = %.1f' % (lo + 1, hi, L, v0, len(terms), umax))
for b in range(len(edges) - 1):
    sel = [(p, m) for p, m in terms if edges[b] <= m*np.log(p)/L < edges[b + 1]]
    if not sel: continue
    # theory: sum c_n^2 alpha_n^2 / 2 (the wave of n has variance alpha_n^2/2 over the zeros, up to the Landau-Gonek term)
    cum_th += sum(cn(p, m, L)**2*al(p, m)**2/2 for p, m in sel)
    # data: remove the fitted waves of this band one by one from the residual (waves are nearly orthogonal)
    for i in range(0, len(sel), 400):
        chunk = sel[i:i + 400]
        fr = np.array([m*np.log(p) for p, m in chunk]); am = np.array([al(p, m) for p, m in chunk])
        W = -np.sin(np.outer(g, fr))*am[None, :]
        coef = (W*res[:, None]).sum(axis=0)/(W*W).sum(axis=0); res = res - W @ coef
    print('  u < %.1f (n <= %7d): explained variance data %.5f | theory %.5f | ratio %.3f | fraction of Var S: data %.1f %%, theory %.1f %%' % (edges[b + 1], int(np.exp(edges[b + 1]*L)), v0 - res.var(), cum_th, (v0 - res.var())/cum_th, 100*(1 - res.var()/v0), 100*cum_th/v0), flush=True)
