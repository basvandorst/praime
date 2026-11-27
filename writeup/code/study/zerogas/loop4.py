import sys, numpy as np, mpmath as mp
from sympy import primerange
from sine import p1
from loop3 import S_mid, OFF
f = int(sys.argv[1]); ymax = int(float(sys.argv[2]))
x_, w_ = np.polynomial.legendre.leggauss(100); sg = (x_ + 1)*3.25; wg = w_*3.25; P1 = np.array([p1(v) for v in sg])
beta = lambda u: float(wg @ (P1*np.cos(np.pi*u*sg)))
T0, x, S = S_mid(f); L = float(mp.log(T0/(2*mp.pi))); S = S - S.mean()
# blocks for the noise: 10 consecutive blocks of midpoints
edges = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8]
bins = {}
ps = np.array(list(primerange(2, ymax + 1)), dtype=np.int64)
u_all = np.log(ps)/L
for b in range(len(edges) - 1):
    sel = ps[(u_all >= edges[b]) & (u_all < edges[b + 1])]
    if len(sel) == 0: continue
    W = np.zeros(len(x))
    for i in range(0, len(sel), 2000):
        chunk = sel[i:i + 2000]
        lp = np.log(chunk.astype(float)); ph0 = np.array([float(mp.fmod(T0*mp.log(int(p)), 2*mp.pi)) for p in chunk])
        W -= (np.sin(ph0[None, :] + x[:, None]*lp[None, :])/np.sqrt(chunk.astype(float))[None, :]).sum(axis=1)/np.pi
    slope = float(S @ W)/float(W @ W)
    blk = np.array_split(np.arange(len(x)), 10); sl = np.array([float(S[i] @ W[i])/float(W[i] @ W[i]) for i in blk])
    um = float(np.mean(np.log(sel.astype(float)))/L)
    print('zeros%d u in [%.2f, %.2f): %6d primes, mean u %.3f: slope %+.3f +- %.3f (blocks) | beta(u) = %+.3f | z = %+.1f' % (f, edges[b], edges[b + 1], len(sel), um, slope, sl.std(ddof=1)/np.sqrt(10), beta(um), (slope - beta(um))/(sl.std(ddof=1)/np.sqrt(10))), flush=True)
