import sys, numpy as np, mpmath as mp
from sympy import primerange
from sine import p1
from loop3 import load, OFF
f = int(sys.argv[1]); ymax = int(float(sys.argv[2]))
T0, n0 = OFF[f]; g = load(f); L = float(mp.log(T0/(2*mp.pi)))
th0 = mp.siegeltheta(T0); d1 = float(mp.log(T0/(2*mp.pi))/2); d2 = float(1/(2*T0)); k0 = mp.floor(th0/mp.pi); frac0 = float(th0/mp.pi - k0)
Sz = ((n0 + 1 - int(k0)) + np.arange(len(g))) - 0.5 - (frac0 + (d1*g + d2*g*g/2)/np.pi) - 1          # S at the zeros (mean of the two sides)
Sz = Sz - Sz.mean(); s = np.diff(g)*L/(2*np.pi)
x_, w_ = np.polynomial.legendre.leggauss(100); sg = (x_ + 1)*3.25; wg = w_*3.25; P1 = np.array([p1(v) for v in sg])
beta = lambda u: float(wg @ (P1*np.cos(np.pi*u*sg)))
edges = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6]
ps = np.array(list(primerange(2, ymax + 1)), dtype=np.int64); u_all = np.log(ps)/L
for b in range(len(edges) - 1):
    sel = ps[(u_all >= edges[b]) & (u_all < edges[b + 1])]
    if len(sel) == 0: continue
    W = np.zeros(len(g))
    for i in range(0, len(sel), 2000):
        chunk = sel[i:i + 2000]; lp = np.log(chunk.astype(float)); ph0 = np.array([float(mp.fmod(T0*mp.log(int(p)), 2*mp.pi)) for p in chunk])
        W -= (np.sin(ph0[None, :] + g[:, None]*lp[None, :])/np.sqrt(chunk.astype(float))[None, :]).sum(axis=1)/np.pi
    blk = np.array_split(np.arange(len(g)), 10); c = float(Sz @ W)/float(W @ W)
    cb = np.array([float(Sz[i] @ W[i])/float(W[i] @ W[i]) for i in blk]); um = float(np.mean(np.log(sel.astype(float)))/L)
    print('zeros%d u in [%.2f, %.2f): mean u %.3f: c(u) at the zeros = %+.3f +- %.3f | c(u) beta(u) = %+.3f | 1 - u = %.3f' % (f, edges[b], edges[b + 1], um, c, cb.std(ddof=1)/np.sqrt(10), c*beta(um), 1 - um), flush=True)
