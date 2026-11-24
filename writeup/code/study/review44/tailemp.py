"""Finite-H long-range tail from EMPIRICAL insertion responses r_u (rough classes mod Q) and HL pair factors S_{>y}(h):
W_tail(H) = (2/((C-1)^2 lam)) (Q/phi(Q)) sum_{h<H} (1-h/H) w2(h) A_W(h mod Q) (S_{>y}(h) - 1),  A_W(h) = (1/phi(Q)) sum_u <r_u, r_{u+h}>."""
import numpy as np, math, json
WIN = [6*10**9 + i*10**8 for i in range(24)]
lam = np.mean([math.log(x + 5e7) for x in WIN])
QR = 510510
okR = np.ones(QR, bool)
for p in (2,3,5,7,11,13,17): okR[::p] = False
RU = np.flatnonzero(okR); NR = len(RU)
def sp(n):
    s = np.ones(n+1, bool); s[:2] = False
    for i in range(2, int(n**0.5)+1):
        if s[i]: s[i*i::i] = False
    return np.flatnonzero(s)
PR = sp(10**7)
def S_tail(hmax, y):
    p = PR[PR > y].astype(float); Ct = math.exp(np.sum(np.log1p(-1/(p-1)**2)) - 1/(1e7*math.log(1e7)))
    f = np.ones(hmax+1)
    for pp in PR[(PR > y) & (PR <= hmax)].tolist(): f[pp::pp] *= (pp-1)/(pp-2)
    return f*Ct
def red(q): return [a for a in range(q) if math.gcd(a, q) == 1]
exec(open('resp.py').read().split('lam = np.mean')[0].split('def red')[1].join(['def red', '']) if False else '')
def within_basis(P, pi):
    C = len(pi); R = np.zeros((C*C, C*C)); V = np.zeros((C*C, C))
    for a in range(C):
        s = np.sqrt(P[a]); R[a*C:(a+1)*C, a*C:(a+1)*C] = np.eye(C) - np.outer(s, s)
        for b in range(C): V[a*C+b, b] = math.sqrt(pi[a]*P[a, b])
    PW = R - R @ V @ np.linalg.pinv(V.T @ R @ V) @ V.T @ R
    w, U = np.linalg.eigh(PW); return U[:, w > 0.5]
def yvec(P, pi):
    C = len(pi); Y = np.zeros((C, C, C*C))
    for a in range(C):
        for b in range(C):
            for b2 in range(C): Y[a, b, a*C+b2] = ((1.0 if b2 == b else 0.0) - P[a, b2]) / math.sqrt(pi[a]*P[a, b2])
    return Y
HS = [10**4, 2*10**4, 5*10**4, 10**5, 2*10**5, 5*10**5, 10**6]
Hmax = max(HS); h = np.arange(1, Hmax+1)
out = {}
for y, Q in ((11, 2310), (13, 30030), (17, 510510)):
    St = S_tail(Hmax, y)
    for q in (3, 5, 7):
        r = red(q); C = len(r); lut = -np.ones(q, np.int64); lut[r] = np.arange(C)
        T = sum(np.load('out/w_%d.npz' % x)['T_%d_1' % q] for x in WIN).astype(float)
        P = T / T.sum(1, keepdims=True); pi = T.sum(1) / T.sum()
        F = yvec(P, pi) @ within_basis(P, pi); D = F.shape[2]
        cnt = sum(np.load('out/w_%d.npz' % x)['rc_comp_%d' % q] + np.load('out/w_%d.npz' % x)['rc_prime_%d' % q] for x in WIN).reshape(NR, C, C).astype(float)
        g = RU % Q; uq, inv = np.unique(g, return_inverse=True)
        agg = np.zeros((len(uq), C, C)); np.add.at(agg, inv, cnt)
        cx = lut[uq % q]
        rq = F[:, cx, :].transpose(1, 0, 2)[:, :, None, :] + F[cx][:, None, :, :] - F[None, :, :, :]
        n = agg.sum((1, 2)); m1 = np.einsum('ulr,ulrd->ud', agg, rq) / n[:, None]
        m2 = np.einsum('ulr,ulrd->u', agg, rq**2) / n; noise = (m2 - (m1**2).sum(1)) / (n - 1)
        A = np.zeros(Q)
        for d in range(D):
            full = np.zeros(Q); full[uq] = m1[:, d]; Fq = np.fft.rfft(full); A += np.fft.irfft(Fq*np.conj(Fq), n=Q)
        A /= len(uq); A[0] -= noise.mean()          # remove sampling noise from the zero lag
        pref = (Q/len(uq)) / ((C-1)**2 * lam)
        base = A[h % Q] * (St[1:Hmax+1] - 1) * (1 - np.exp(-h/lam))
        Wt = {H: float(2*pref*np.sum((1 - h[:H-1]/H) * base[:H-1])) for H in HS}
        out['%d_%d' % (q, y)] = dict(V_per_dir=float(A[0])/(C-1)**2, W_tail=Wt)
        print('y=%2d q=%d V/dir %.4f  W_tail: %s   slopes/decade 1e4->1e5 %+.4f 1e5->1e6 %+.4f' % (y, q, A[0]/(C-1)**2, ' '.join('%+.4f' % Wt[H] for H in HS), Wt[10**5]-Wt[10**4], Wt[10**6]-Wt[10**5]), flush=True)
json.dump(out, open('tailemp.json', 'w'), indent=1)
