import sys, numpy as np
from sympy import primerange
Z1 = '/home/bas/praime/writeup/data/riemann/odlyzko_zeros1.txt'
lo, hi = int(sys.argv[1]), int(sys.argv[2]); umax = float(sys.argv[3]) if len(sys.argv) > 3 else 1.2
g_all = np.loadtxt(Z1); g = g_all[lo:hi]                      # zeros number lo+1 .. hi
theta = lambda t: t/2*np.log(t/(2*np.pi)) - t/2 - np.pi/8 + 1/(48*t) + 7/(5760*t**3)
n = np.arange(lo + 1, hi + 1)
Sz = n - 0.5 - (theta(g)/np.pi + 1)                           # S at the zeros (mean of both sides); N(t) = theta/pi + 1 + S
mid = (g[:-1] + g[1:])/2; Sm = n[:-1] - (theta(mid)/np.pi + 1)
print('zeros %d..%d: t from %.0f to %.0f; mean S at zeros %.4f, at midpoints %.4f; var %.4f, %.4f' % (lo + 1, hi, g[0], g[-1], Sz.mean(), Sm.mean(), Sz.var(), Sm.var()))
Lz = np.log(g/(2*np.pi)); Lm = np.log(mid/(2*np.pi)); Lc = float(Lz.mean())
edges = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]
nmax = int(np.exp(umax*Lz.max())) + 1
# prime powers: frequency ln(p^m), amplitude 1/(m p^(m/2))
terms = []
for p in primerange(2, nmax + 1):
    m = 1
    while p**m <= nmax: terms.append((m*np.log(p), 1/(m*p**(m/2)))); m += 1
fr = np.array([t[0] for t in terms]); am = np.array([t[1] for t in terms]); u = fr/Lc
cue = np.load('cue_N400.npy'); cue_edges = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2]
def cue_val(row, a, b):
    idx = [i for i in range(len(cue_edges) - 1) if cue_edges[i] >= a - 1e-9 and cue_edges[i + 1] <= b + 1e-9]
    return float(np.mean(cue[row][idx])) if idx else float('nan')
Szc = Sz - Sz.mean(); Smc = Sm - Sm.mean()
for b in range(len(edges) - 1):
    sel = (u >= edges[b]) & (u < edges[b + 1])
    if not sel.any(): continue
    Wz = np.zeros(len(g)); Wm = np.zeros(len(mid))
    idx = np.nonzero(sel)[0]
    for i in range(0, len(idx), 500):
        j = idx[i:i + 500]
        Wz -= (np.sin(np.outer(g, fr[j]))*am[j][None, :]).sum(axis=1)/np.pi
        Wm -= (np.sin(np.outer(mid, fr[j]))*am[j][None, :]).sum(axis=1)/np.pi
    blk = np.array_split(np.arange(len(mid)), 10)
    c = float(Szc @ Wz)/float(Wz @ Wz); bb = float(Smc @ Wm)/float(Wm @ Wm)
    ce = np.array([float(Szc[i] @ Wz[i])/float(Wz[i] @ Wz[i]) for i in blk]).std(ddof=1)/np.sqrt(10)
    be = np.array([float(Smc[i] @ Wm[i])/float(Wm[i] @ Wm[i]) for i in blk]).std(ddof=1)/np.sqrt(10)
    print('u in [%.1f, %.1f): %6d prime powers: c = %+.3f +- %.3f (matrices %+.3f) | b = %+.3f +- %.3f (matrices %+.3f)' % (edges[b], edges[b + 1], sel.sum(), c, ce, cue_val(0, edges[b], edges[b + 1]), bb, be, cue_val(2, edges[b], edges[b + 1])), flush=True)
