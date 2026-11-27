import sys, numpy as np
from sympy import primerange
from loop11 import cn, al
lo, hi, umax = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
g = np.loadtxt('/home/bas/praime/writeup/data/riemann/odlyzko_zeros1.txt')[lo:hi]
theta = lambda t: t/2*np.log(t/(2*np.pi)) - t/2 - np.pi/8 + 1/(48*t) + 7/(5760*t**3)
Sz = np.arange(lo + 1, hi + 1) - 0.5 - (theta(g)/np.pi + 1); Sz = Sz - Sz.mean(); Lz = np.log(g/(2*np.pi)); L = float(Lz.mean())
nmax = int(np.exp(umax*L)) + 1; terms = []
for p in primerange(2, nmax + 1):
    m = 1
    while p**m <= nmax: terms.append((p, m)); m += 1
terms.sort(key=lambda pm: pm[1]*np.log(pm[0]))
edges = np.arange(0, umax + 1e-9, 0.1); pred = np.zeros(len(g)); cum = 0.0; v0 = Sz.var()
print('zeros %d..%d: L = %.3f, Var S at the zeros = %.5f' % (lo + 1, hi, L, v0))
for b in range(len(edges) - 1):
    sel = [(p, m) for p, m in terms if edges[b] <= m*np.log(p)/L < edges[b + 1]]
    if not sel: continue
    for i in range(0, len(sel), 400):
        chunk = sel[i:i + 400]
        fr = np.array([m*np.log(p) for p, m in chunk]); am = np.array([al(p, m) for p, m in chunk]); c = np.array([cn(p, m, L) for p, m in chunk])
        pred -= np.sin(np.outer(g, fr)) @ (am*c)
        cum += float(np.sum(c*c*am*am/2*(1 + 2*np.pi*fr*np.array([al(p, 2*m) for p, m in chunk])/L)))
    res = Sz - pred
    print('  n <= %7d (u < %.1f): variance removed %.5f | predicted %.5f | ratio %.3f | residual %.1f %% of Var S, rms error of a zero %.3f spacings' % (int(np.exp(edges[b + 1]*L)), edges[b + 1], v0 - res.var(), cum, (v0 - res.var())/cum, 100*res.var()/v0, np.sqrt(res.var())), flush=True)
