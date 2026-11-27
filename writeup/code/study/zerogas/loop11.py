import sys, numpy as np, mpmath as mp
from loop9 import cth
mp.mp.dps = 50
PP = [(2,1),(3,1),(2,2),(5,1),(7,1),(2,3),(3,2),(11,1),(13,1),(2,4),(17,1),(19,1),(23,1),(5,2),(3,3),(29,1),(31,1),(2,5)]
al = lambda p, m: 1/(np.pi*m*p**(m/2))
def cn(p, m, L):
    a = al(p, m); f = m*np.log(p)
    prod = sum(al(p, i)*al(p, m - i) for i in range(1, m)); diff = sum(al(p, m + j)*al(p, j) for j in range(1, 60))
    Delta = -(a*f/8)*(prod - 2*diff)
    return (float(cth(np.array([f/L]))[0]) + 4*np.pi*Delta/(L*a*a))/(1 + 2*np.pi*f*al(p, 2*m)/L)
def run(g, Sz, L, phase0, label):
    Szc = Sz - Sz.mean(); blk = np.array_split(np.arange(len(g)), 10)
    for (p, m), ph in zip(PP, phase0):
        f = m*np.log(p); W = -al(p, m)*np.sin(ph + g*f)
        c = float(Szc @ W)/float(W @ W); e = np.array([float(Szc[i] @ W[i])/float(W[i] @ W[i]) for i in blk]).std(ddof=1)/np.sqrt(10)
        pred = cn(p, m, L); mat = float(cth(np.array([f/L]))[0])
        print('%s n = %d^%d = %3d (u = %.3f): c = %+.4f +- %.4f | matrices %+.4f | with the arithmetic correction %+.4f | z = %+.1f' % (label, p, m, p**m, f/L, c, e, mat, pred, (c - pred)/e), flush=True)
if __name__ == '__main__':
    which = sys.argv[1]
    if which in ('3', '4', '5'):
        from loop3 import load, OFF
        f_ = int(which); T0, n0 = OFF[f_]; g = load(f_); L = float(mp.log(T0/(2*mp.pi)))
        th0 = mp.siegeltheta(T0); d1 = float(mp.log(T0/(2*mp.pi))/2); d2 = float(1/(2*T0)); k0 = mp.floor(th0/mp.pi); frac0 = float(th0/mp.pi - k0)
        Sz = ((n0 + 1 - int(k0)) + np.arange(len(g))) - 0.5 - (frac0 + (d1*g + d2*g*g/2)/np.pi) - 1
        run(g, Sz, L, [float(mp.fmod(T0*m*mp.log(p), 2*mp.pi)) for p, m in PP], 'zeros%d' % f_)
    else:
        lo, hi = int(sys.argv[2]), int(sys.argv[3])
        g = np.loadtxt('/home/bas/praime/writeup/data/riemann/odlyzko_zeros1.txt')[lo:hi]
        theta = lambda t: t/2*np.log(t/(2*np.pi)) - t/2 - np.pi/8 + 1/(48*t) + 7/(5760*t**3)
        Sz = np.arange(lo + 1, hi + 1) - 0.5 - (theta(g)/np.pi + 1); L = float(np.log(g/(2*np.pi)).mean())
        run(g, Sz, L, [0.0]*len(PP), 'zeros %d..%d' % (lo + 1, hi))
