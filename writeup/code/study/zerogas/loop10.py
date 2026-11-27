import sys, numpy as np, mpmath as mp
from sympy import primerange
from loop9 import cth
mp.mp.dps = 50
def bands(g, Sz, L, terms_fr, terms_am, phase0, edges, label):
    """g: heights relative to an offset; phase0: (offset * frequency) mod 2 pi for each term."""
    u = terms_fr/L; Szc = Sz - Sz.mean(); blk = np.array_split(np.arange(len(g)), 10); zs = []
    for b in range(len(edges) - 1):
        idx = np.nonzero((u >= edges[b]) & (u < edges[b + 1]))[0]
        if len(idx) == 0: continue
        W = np.zeros(len(g))
        for i in range(0, len(idx), 1000):
            j = idx[i:i + 1000]; W -= (np.sin(phase0[j][None, :] + np.outer(g, terms_fr[j]))*terms_am[j][None, :]).sum(axis=1)/np.pi
        c = float(Szc @ W)/float(W @ W); e = np.array([float(Szc[i] @ W[i])/float(W[i] @ W[i]) for i in blk]).std(ddof=1)/np.sqrt(10)
        wts = terms_am[idx]**2; th = float((cth(u[idx])*wts).sum()/wts.sum()); zs.append((c - th)/e)
        print('%s u in [%.2f, %.2f): %6d prime powers, c = %+.4f +- %.4f | formula %+.4f | z = %+.1f' % (label, edges[b], edges[b + 1], len(idx), c, e, th, zs[-1]), flush=True)
    return zs
def terms(nmax):
    fr, am, pw = [], [], []
    for p in primerange(2, nmax + 1):
        m = 1
        while p**m <= nmax: fr.append(m*np.log(p)); am.append(1/(m*p**(m/2))); pw.append((p, m)); m += 1
    return np.array(fr), np.array(am), pw
which = sys.argv[1]
if which in ('3', '4', '5'):
    from loop3 import load, OFF
    f = int(which); T0, n0 = OFF[f]; g = load(f); L = float(mp.log(T0/(2*mp.pi)))
    th0 = mp.siegeltheta(T0); d1 = float(mp.log(T0/(2*mp.pi))/2); d2 = float(1/(2*T0)); k0 = mp.floor(th0/mp.pi); frac0 = float(th0/mp.pi - k0)
    Sz = ((n0 + 1 - int(k0)) + np.arange(len(g))) - 0.5 - (frac0 + (d1*g + d2*g*g/2)/np.pi) - 1
    fr, am, pw = terms(int(float(sys.argv[2])))
    ph0 = np.array([float(mp.fmod(T0*m*mp.log(p), 2*mp.pi)) for p, m in pw])
    bands(g, Sz, L, fr, am, ph0, np.arange(0, 0.6001, 0.05), 'zeros%d' % f)
else:
    lo, hi = int(sys.argv[2]), int(sys.argv[3]); umax = float(sys.argv[4])
    g = np.loadtxt('/home/bas/praime/writeup/data/riemann/odlyzko_zeros1.txt')[lo:hi]
    theta = lambda t: t/2*np.log(t/(2*np.pi)) - t/2 - np.pi/8 + 1/(48*t) + 7/(5760*t**3)
    Sz = np.arange(lo + 1, hi + 1) - 0.5 - (theta(g)/np.pi + 1); L = float(np.log(g/(2*np.pi)).mean())
    fr, am, pw = terms(int(np.exp(umax*np.log(g.max()/(2*np.pi)))) + 1)
    bands(g, Sz, L, fr, am, np.zeros(len(fr)), np.arange(0, umax + 1e-9, 0.1), 'zeros %d..%d' % (lo + 1, hi))
