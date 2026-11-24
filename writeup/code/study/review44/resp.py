"""Empirical insertion responses r_u per rough class (counterfactual: neighbours of each rough site), V = mean_u |Pi_W r_u|^2."""
import numpy as np, math, json, sys
WIN = [6*10**9 + i*10**8 for i in range(24)]
QR = 510510
ok = np.ones(QR, bool)
for p in (2,3,5,7,11,13,17): ok[::p] = False
RU = np.flatnonzero(ok); NR = len(RU)
def red(q): return [a for a in range(q) if math.gcd(a, q) == 1]
def within_basis(P, pi):
    C = len(pi); R = np.zeros((C*C, C*C)); V = np.zeros((C*C, C))
    for a in range(C):
        s = np.sqrt(P[a]); R[a*C:(a+1)*C, a*C:(a+1)*C] = np.eye(C) - np.outer(s, s)
        for b in range(C): V[a*C+b, b] = math.sqrt(pi[a]*P[a, b])
    PW = R - R @ V @ np.linalg.pinv(V.T @ R @ V) @ V.T @ R
    w, U = np.linalg.eigh(PW); U = U[:, w > 0.5]
    assert U.shape[1] == (C-1)**2, U.shape
    return U
def yvec(P, pi):
    C = len(pi); Y = np.zeros((C, C, C*C))
    for a in range(C):
        for b in range(C):
            for b2 in range(C): Y[a, b, a*C+b2] = ((1.0 if b2 == b else 0.0) - P[a, b2]) / math.sqrt(pi[a]*P[a, b2])
    return Y
lam = np.mean([math.log(x + 5e7) for x in WIN])
res = {}
for q in (3, 5, 7):
    r = red(q); C = len(r); lut = -np.ones(q, np.int64); lut[r] = np.arange(C)
    T = sum(np.load('out/w_%d.npz' % x)['T_%d_1' % q] for x in WIN).astype(float)
    P = T / T.sum(1, keepdims=True); pi = T.sum(1) / T.sum()
    U = within_basis(P, pi); F = yvec(P, pi) @ U            # (C, C, D) within coords of y(a,b)
    D = F.shape[2]
    cc = sum(np.load('out/w_%d.npz' % x)['rc_comp_%d' % q] for x in WIN).reshape(NR, C, C).astype(float)
    cp = sum(np.load('out/w_%d.npz' % x)['rc_prime_%d' % q] for x in WIN).reshape(NR, C, C).astype(float)
    cx = lut[RU % q]
    # response per (u, cL, cR): F[cL,cx] + F[cx,cR] - F[cL,cR]
    resp = F[:, cx, :].transpose(1, 0, 2)[:, :, None, :] + F[cx][:, None, :, :] - F[None, :, :, :]   # (NR, C, C, D)
    out = {}
    for name, cnt in (('all', cc + cp), ('comp', cc), ('prime', cp)):
        for Qp in (210, 2310, 30030, 510510):
            if Qp % q: continue
            g = RU % Qp
            # aggregate counts by class mod Qp
            uq, inv = np.unique(g, return_inverse=True)
            agg = np.zeros((len(uq), C, C)); np.add.at(agg, inv, cnt)
            cxq = lut[uq % q]
            rq = F[:, cxq, :].transpose(1, 0, 2)[:, :, None, :] + F[cxq][:, None, :, :] - F[None, :, :, :]
            n = agg.sum((1, 2))
            m1 = np.einsum('ulr,ulrd->ud', agg, rq) / n[:, None]
            m2 = np.einsum('ulr,ulrd->u', agg, rq ** 2) / n
            var = m2 - (m1 ** 2).sum(1)
            Vraw = float(((m1 ** 2).sum(1)).mean()); bias = float((var / (n - 1)).mean())
            out['%s_%d' % (name, Qp)] = dict(V_per_dir=(Vraw - bias) / (C - 1) ** 2, bias_per_dir=bias / (C - 1) ** 2, ncls=len(uq), nmin=float(n.min()),
                                              mean_resp_norm2=float((np.einsum('ulr,ulrd->d', agg, rq) / n.sum()) @ (np.einsum('ulr,ulrd->d', agg, rq) / n.sum())) / (C - 1) ** 2)
            print('q=%d %-5s Q=%6d: V/dir = %.4f (noise bias %.5f, %d classes, min n %d); slope -V/lam per unit lnH = %.5f, per decade %.4f' % (q, name, Qp, out['%s_%d' % (name, Qp)]['V_per_dir'], bias / (C-1)**2, len(uq), n.min(), -out['%s_%d' % (name, Qp)]['V_per_dir'] / lam, -out['%s_%d' % (name, Qp)]['V_per_dir'] / lam * math.log(10)), flush=True)
    res[q] = out
    np.save('rmean_%d.npy' % q, None) if False else None
json.dump(dict(lam=lam, res=res), open('resp.json', 'w'), indent=1)
