import sys, numpy as np
from sympy import primerange
from loop9 import cth
from loop11 import al
lo, hi, umax = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
g = np.loadtxt('/home/bas/praime/writeup/data/riemann/odlyzko_zeros1.txt')[lo:hi]
theta = lambda t: t/2*np.log(t/(2*np.pi)) - t/2 - np.pi/8 + 1/(48*t) + 7/(5760*t**3)
Sz = np.arange(lo + 1, hi + 1) - 0.5 - (theta(g)/np.pi + 1); Sz = Sz - Sz.mean(); Lz = np.log(g/(2*np.pi)); L = float(Lz.mean())
nmax = int(np.exp(umax*Lz.max())) + 1; terms = []
for p in primerange(2, nmax + 1):
    m = 1
    while p**m <= nmax: terms.append((p, m)); m += 1
terms.sort(key=lambda pm: pm[1]*np.log(pm[0]))
def parts(p, m):
    a = al(p, m); f = m*np.log(p)
    prod = sum(al(p, i)*al(p, m - i) for i in range(1, m)); diff = sum(al(p, m + j)*al(p, j) for j in range(1, 60))
    return a, f, 4*np.pi*(-(a*f/8)*(prod - 2*diff))/(a*a), 2*np.pi*f*al(p, 2*m)       # amplitude, frequency, numerator term * L, denominator term * L
pred = np.zeros(len(g)); v0 = Sz.var(); edges = np.arange(0, umax + 1e-9, 0.1)
print('zeros %d..%d: Var S at the zeros = %.5f' % (lo + 1, hi, v0))
for b in range(len(edges) - 1):
    sel = [(p, m) for p, m in terms if edges[b] <= m*np.log(p)/L < edges[b + 1]]
    for i in range(0, len(sel), 200):
        P = [parts(p, m) for p, m in sel[i:i + 200]]
        a = np.array([x[0] for x in P]); f = np.array([x[1] for x in P]); nu = np.array([x[2] for x in P]); de = np.array([x[3] for x in P])
        U = f[None, :]/Lz[:, None]
        C = (cth(U.ravel()).reshape(U.shape)*(U < 2) + nu[None, :]/Lz[:, None])/(1 + de[None, :]/Lz[:, None])
        C = np.where(U < umax, C, 0.0)                              # a zero uses a prime power only while its own u is below umax
        pred -= (np.sin(np.outer(g, f))*C) @ a
    res = Sz - pred
    print('  u < %.1f: residual %.2f %% of Var S, rms error of a zero %.4f spacings' % (edges[b + 1], 100*res.var()/v0, np.sqrt(res.var())), flush=True)
