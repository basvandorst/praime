import sys, numpy as np
from sine import p1
rng = np.random.default_rng(16); N, M = 400, 300
ks = np.arange(1, int(0.6*N) + 1); num = np.zeros((M, len(ks))); den = np.zeros((M, len(ks)))
for it in range(M):
    Z = (rng.standard_normal((N, N)) + 1j*rng.standard_normal((N, N)))/np.sqrt(2); Q, R = np.linalg.qr(Z)
    x = np.sort(np.angle(np.linalg.eigvals(Q*(np.diag(R)/np.abs(np.diag(R)))[None, :])) % (2*np.pi))*N/(2*np.pi)
    Sz = np.arange(N) - x; dS = np.roll(Sz, -1) - Sz; dS[-1] += 0.0          # S(x_{n+1}) - S(x_n) = 1 - s_n (circular: last one wraps)
    dS[-1] = 1 - (x[0] + N - x[-1])
    E = np.exp(-2j*np.pi*np.outer(ks, x)/N); coef = E.sum(axis=1)/(2j*np.pi*ks); W = 2*np.real(coef[:, None]*np.conj(E))
    dW = np.roll(W, -1, axis=1) - W
    num[it] = dW @ dS; den[it] = (dW*dW).sum(axis=1)
u = ks/N; cp = num.sum(0)/den.sum(0)
# smooth c'(u) in bands of 0.02 and report
edges = np.arange(0, 0.6001, 0.02); mid = []; val = []; err = []
grp = np.array_split(np.arange(M), 10)
for i in range(len(edges) - 1):
    sel = (u >= edges[i]) & (u < edges[i + 1])
    if not sel.any(): continue
    mid.append(u[sel].mean()); val.append(num[:, sel].sum()/den[:, sel].sum())
    err.append(np.array([num[g][:, sel].sum()/den[g][:, sel].sum() for g in grp]).std(ddof=1)/np.sqrt(10))
mid, val, err = map(np.array, (mid, val, err))
print("c'(u) from matrices:", ' '.join('%.2f:%.4f' % (a, b) for a, b in zip(mid[::2], val[::2])))
cfun = lambda uu: np.interp(uu, mid, val)
# spacing law of the sine process, for the weight of each power: E[(W(x+s) - W(x))^2] = alpha^2 (1 - E cos(2 pi u s))
x_, w_ = np.polynomial.legendre.leggauss(100); sg = (x_ + 1)*3.25; wg = w_*3.25; P1 = np.array([p1(v) for v in sg])
for L in (10.23, 10.5, 11.13):
    out = []
    for p in (2, 3, 5):
        tot = 0.0; acc = 0.0
        for k in range(1, 12):
            uu = k*np.log(p)/L
            if uu >= 0.58: break
            v = (1/(k*k*p**k))*(1 - float(wg @ (P1*np.cos(2*np.pi*uu*sg))))
            tot += v; acc += v*(2*cfun(uu) - 1)
        out.append(acc/tot)
    d = np.array([0.0126, 0.0144, 0.0136]); amp = float((np.array(out)*d*d).sum()/(d*d).sum())
    print("L = %.2f: amplitude against the closure for p = 2, 3, 5: %.4f %.4f %.4f; weighted as in round 147: %.4f (measured 0.940 +- 0.011)" % (L, out[0], out[1], out[2], amp))
