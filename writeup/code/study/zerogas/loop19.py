import sys, numpy as np
from sympy import primerange, factorint
from loop9 import cth
from loop11 import al
lo, hi, ucut, xmax = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
MODE = sys.argv[5] if len(sys.argv) > 5 else 'law'
g = np.loadtxt('/home/bas/praime/writeup/data/riemann/odlyzko_zeros1.txt')[lo:hi]
theta = lambda t: t/2*np.log(t/(2*np.pi)) - t/2 - np.pi/8 + 1/(48*t) + 7/(5760*t**3)
idx = np.arange(lo + 1, hi + 1).astype(float)
# prime powers used: u = ln n / L(t) < ucut at the START of the window (so every zero has them below ucut)
Lmin = np.log(g[0]/(2*np.pi)); nmax = int(np.exp(ucut*Lmin)); terms = []
for p in primerange(2, nmax + 1):
    m = 1
    while p**m <= nmax:
        a = al(p, m); f = m*np.log(p)
        prod = sum(al(p, i)*al(p, m - i) for i in range(1, m)); diff = sum(al(p, m + j)*al(p, j) for j in range(1, 60))
        terms.append((a, f, 4*np.pi*(-(a*f/8)*(prod - 2*diff))/(a*a), 2*np.pi*f*al(p, 2*m))); m += 1
A = np.array([t[0] for t in terms]); F = np.array([t[1] for t in terms]); NU = np.array([t[2] for t in terms]); DE = np.array([t[3] for t in terms])
def Spred(y):
    out = np.zeros(len(y))
    for i in range(0, len(y), 2000):
        yy = y[i:i + 2000]; Ly = np.log(yy/(2*np.pi))[:, None]; U = F[None, :]/Ly
        C = (cth(U.ravel()).reshape(U.shape) + NU[None, :]/Ly)/(1 + DE[None, :]/Ly)
        if MODE == 'plain': C = np.ones_like(U)
        if MODE == 'matrix': C = cth(U.ravel()).reshape(U.shape)
        out[i:i + 2000] = -((np.sin(np.outer(yy, F))*C) @ A)
    return out
def solve_smooth(v, y):                       # theta(y)/pi + 1 = v
    for _ in range(8): y = y - (theta(y)/np.pi + 1 - v)/(np.log(y/(2*np.pi))/(2*np.pi))
    return y
y = solve_smooth(idx - 0.5, 2*np.pi*idx/np.log(idx))          # zeros with no prime at all (start)
y0 = y.copy()
for it in range(40):
    ynew = solve_smooth(idx - 0.5 - Spred(y), y); y = 0.5*y + 0.5*ynew
sp = 2*np.pi/np.log(g/(2*np.pi))
print('zeros %d..%d, t = %.0f..%.0f; prime powers n <= %d used (%d terms)' % (lo + 1, hi, g[0], g[-1], nmax, len(terms)))
print('rms error of a zero in mean spacings: no primes %.4f | rebuilt %.4f' % (np.sqrt(np.mean(((y0 - g)/sp)**2)), np.sqrt(np.mean(((y - g)/sp)**2))))
if xmax <= nmax: sys.exit(0)
xs = np.arange(nmax + 1, xmax + 1); lx = np.log(xs.astype(float))
Ft = np.zeros(len(xs)); Fr = np.zeros(len(xs)); F0 = np.zeros(len(xs))
for i in range(0, len(xs), 400):
    l = lx[i:i + 400]
    Ft[i:i + 400] = np.cos(np.outer(l, g)).sum(axis=1); Fr[i:i + 400] = np.cos(np.outer(l, y)).sum(axis=1); F0[i:i + 400] = np.cos(np.outer(l, y0)).sum(axis=1)
lam = np.zeros(len(xs))
for j, x in enumerate(xs):
    fa = factorint(int(x))
    if len(fa) == 1: lam[j] = np.log(list(fa)[0])
target = -lam/np.sqrt(xs)
tm = np.exp(np.log(g/(2*np.pi)).mean())
cor = lambda a, b: float(np.corrcoef(a, b)[0, 1])
for a, b in ((nmax, tm), (tm, 1.5*tm), (1.5*tm, 2*tm), (2*tm, xmax)):
    sel = (xs > a) & (xs <= b)
    if sel.sum() < 50: continue
    print('x in (%.0f, %.0f] (u %.2f..%.2f, %d integers, %d prime powers): corr(rebuilt, true) %.3f | corr with -Lambda/sqrt(x): true zeros %.3f, rebuilt %.3f, no primes %.3f | line amplitude rebuilt/true %.3f' % (
        a, b, np.log(a)/np.log(tm), np.log(b)/np.log(tm), sel.sum(), (lam[sel] > 0).sum(), cor(Fr[sel], Ft[sel]), cor(Ft[sel], target[sel]), cor(Fr[sel], target[sel]), cor(F0[sel], target[sel]),
        float(Fr[sel] @ target[sel])/float(Ft[sel] @ target[sel])), flush=True)
np.save('loop19_spectra.npy', np.array([xs, Ft, Fr, F0, lam]))
