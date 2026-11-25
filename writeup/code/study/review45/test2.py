import numpy as np, time, sys
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
from mylaw import *
nd = nodes(); print('nodes', len(nd[2]))
ls = np.array([1, 2.5, 4, 5.5, 8.0])
D = 1004; bad = badset(D); kp2 = kphi2(D); L = 10.5
for P, tail in ((3000, False), (3000, True), (30000, True)):
    t0 = time.time()
    import mylaw
    f = T3_P(nd[2], nd[0][nd[4]], L, bad, kp2, P=P) if False else None
    S, WS, T, WT, IDX = nd
    fv = np.zeros(len(T))
    c = 2*np.pi/L
    # inline to pass tail flag
    u = np.zeros_like(T); v = c*T; w = c*S[IDX]
    tot = np.zeros(T.shape, complex)
    for (z, x, yv) in ((w, u, v), (v, w, u), (u, w, v)):
        X = x - z; Yv = yv - z
        p1X = P1(X, L, bad, kp2); p1Y = P1(Yv, L, bad, kp2)
        P3a = B1(X, Yv, bad, P, tail) + Llogd(Yv, bad) - Llogd(Yv - X, bad)
        P3b = B1(Yv, X, bad, P, tail) + Llogd(X, bad) - Llogd(X - Yv, bad)
        tot = tot - p1X*P3a - p1Y*P3b
    fv = 2*tot.real/L**3
    F = 2*np.bincount(IDX, weights=fv*WT, minlength=len(S))
    k = np.array([6*np.sum(WS*(l - S)*F*(S < l)) for l in ls])
    print(P, tail, np.round(k, 8), 'time', round(time.time()-t0, 1))
# finer grid check
nd2 = nodes(piece=0.125, ng=8)
t0 = time.time()
print('fine grid', np.round(kf3_P(ls, L, bad, kp2, P=3000, nd=nd2), 8), round(time.time()-t0,1))
