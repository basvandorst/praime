import numpy as np, json, sys, time, os
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
from mylaw import *
SCR = __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/'
ls = np.arange(0.5, 8.01, 0.5)
ND = nodes()
xg, wg = np.polynomial.legendre.leggauss(16)
R = np.concatenate([a + 0.25*(xg + 1)/2 for a in np.arange(0, 8, 0.25)])
WR = np.concatenate([0.25*wg/2 for a in np.arange(0, 8, 0.25)])

def kQ_fast(L, bad):
    pr = PR_ALL[~np.isin(PR_ALL, list(bad))]
    out = np.zeros(len(ls))
    for N in range(2, 61):
        sel = pr[pr**(-float(N)) > 1e-19]
        if len(sel) == 0: break
        th = np.pi*ls[:, None]*np.log(sel)[None, :]/L
        w = sel**(-float(N))
        for m in range(1, N):
            n = N - m
            out += (np.sin(m*th)*np.sin(n*th)*np.sin(N*th)/(m*n*N)*w).sum(1)
    return -(6/np.pi**3)*out

def at_L(L, bad, kp2):
    kp = kf3_P(ls, L, bad, kp2, P=3000, nd=ND)
    kq = kQ_fast(L, bad)
    t2 = T2(R, L, bad, kp2)
    s2 = np.array([l + 2*np.sum(WR*(l - R)*t2*(R < l)) for l in ls])
    return kp, kq, s2

def law(D, nh=6, lo=100., hi=600.):
    bad = badset(D); kp2 = kphi2(D)
    x, w = np.polynomial.legendre.leggauss(nh)
    t = lo + (hi - lo)*(x + 1)/2
    Ls = np.log(abs(D)*t/(2*np.pi))
    W = w*Ls; W = W/W.sum()
    acc = np.zeros((3, len(ls)))
    for Li, Wi in zip(Ls, W):
        acc += Wi*np.array(at_L(Li, bad, kp2))
    kp, kq, s2 = acc
    return dict(D=D, meanL=float((W*Ls).sum()), kfP=kp.tolist(), kQ=kq.tolist(), sigma2=s2.tolist(),
                kf3=(kp + kq).tolist(), full=(kp + kq + 3*s2 - 2*ls).tolist())

if __name__ == '__main__':
    Ds = [int(a) for a in sys.argv[1].split(',')]
    nh = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    fn = SCR + 'mylaw_out.json'
    res = json.load(open(fn)) if os.path.exists(fn) else {}
    for D in Ds:
        t0 = time.time()
        r = law(D, nh)
        res[f'{D}_{nh}'] = r
        json.dump(res, open(fn, 'w'))
        print(D, nh, round(time.time() - t0, 1), np.round(r['full'], 5).tolist(), flush=True)
