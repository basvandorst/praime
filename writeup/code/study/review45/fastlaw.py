import numpy as np, json, sys, time, os
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
from mylaw import *
SCR = __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/'
ls = np.arange(0.5, 8.01, 0.5)
NDF = nodes(piece=0.5, ng=8)
xg, wg = np.polynomial.legendre.leggauss(12)
R = np.concatenate([a + 0.5*(xg + 1)/2 for a in np.arange(0, 8, 0.5)]); WR = np.concatenate([0.5*wg/2 for a in np.arange(0, 8, 0.5)])

def kQ_parts(L, bad):
    """closed-form diagonal: (all m,n over p not | k), (m=n=1 only over p not | k), (all m,n over all primes)"""
    out = np.zeros((3, len(ls)))
    for N in range(2, 61):
        sel = PR_ALL[PR_ALL**(-float(N)) > 1e-19]
        if len(sel) == 0: break
        good = ~np.isin(sel, list(bad))
        th = np.pi*ls[:, None]*np.log(sel)[None, :]/L; w = sel**(-float(N))
        for m in range(1, N):
            n = N - m
            term = np.sin(m*th)*np.sin(n*th)*np.sin(N*th)/(m*n*N)*w
            out[0] += term[:, good].sum(1); out[2] += term.sum(1)
            if N == 2: out[1] += term[:, good].sum(1)
    return -(6/np.pi**3)*out

def at_L(L, bad, kp2, P=1000):
    kp = kf3_P(ls, L, bad, kp2, P=P, nd=NDF)
    y = 2*np.pi*R/L
    P2 = Llogdd(y, bad) - Bsum(y, bad, P=100000)
    t2 = 2*(P1(y, L, bad, kp2) + P2).real/L**2
    s2 = np.array([l + 2*np.sum(WR*(l - R)*t2*(R < l)) for l in ls])
    kq = kQ_parts(L, bad)
    return np.vstack([kp, s2, kq])     # rows: kfP, sigma2, kQ, kQ11, kQall

def band(D, lo, hi, nh=3, shift=0.0):
    bad = badset(D); kp2 = kphi2(D)
    x, w = np.polynomial.legendre.leggauss(nh)
    t = lo + (hi - lo)*(x + 1)/2; Ls = np.log(abs(D)*t/(2*np.pi)); W = w*Ls; W = W/W.sum()
    acc = 0
    for Li, Wi in zip(Ls, W): acc = acc + Wi*at_L(Li + shift, bad, kp2)
    return acc

if __name__ == '__main__':
    if sys.argv[1] == 'test':
        for D in (1004, -1155):
            t0 = time.time(); b = band(D, 100, 600, nh=6); dt = time.time() - t0
            full = b[0] + b[2] + 3*b[1] - 2*ls
            my = json.load(open(SCR + 'mylaw_out.json'))[f'{D}_6']
            print(D, 'fast vs slow max diff full %.1e' % np.max(np.abs(full - np.array(my['full']))), 'time per node %.2f s' % (dt/6))
        sys.exit()
    P = json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/triple/predictions_kappa3.json'))
    Ds = P['Fn']['D']; i0, i1 = int(sys.argv[1]), int(sys.argv[2])
    fn = SCR + 'fastlaw_out.json'
    res = json.load(open(fn)) if os.path.exists(fn) else {}
    t0 = time.time()
    for D in Ds[i0:i1]:
        if str(D) in res: continue
        r = {'lo': band(D, 100, 350, nh=2).tolist(), 'hi': band(D, 350, 600, nh=2).tolist(),
             'Lm': band(D, 100, 600, nh=2, shift=-0.3).tolist(), 'Lp': band(D, 100, 600, nh=2, shift=+0.3).tolist()}
        res[str(D)] = r
        json.dump(res, open(fn, 'w'))
    print('done', i0, i1, round(time.time() - t0, 1))
