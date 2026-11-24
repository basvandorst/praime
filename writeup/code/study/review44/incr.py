"""Explained response variance V(Q) = mean over rough u mod Q of |Pi_W r_u|^2 (noise-corrected), Q = 210 and 210p;
increment inc(p) = V(210p) - V(210): how much residue information mod p adds."""
import numpy as np, math, json
WIN = [6*10**9 + i*10**8 for i in range(24)]
lam = np.mean([math.log(x+5e7) for x in WIN])
PL = [11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 151, 199, 251, 307, 401, 503]
def red(q): return [a for a in range(q) if math.gcd(a, q) == 1]
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
out = {}
for q in (3, 5, 7):
    r = red(q); C = len(r); lut = -np.ones(q, np.int64); lut[r] = np.arange(C)
    T = sum(np.load('out/w_%d.npz' % x)['T_%d_1' % q] for x in WIN).astype(float)
    P = T/T.sum(1, keepdims=True); pi = T.sum(1)/T.sum(); F = yvec(P, pi) @ within_basis(P, pi)
    def V_of(cnt, M, extra_p=None):
        cnt = cnt.reshape(M, C, C).astype(float); u = np.arange(M)
        keep = np.gcd(u, 210) == 1
        if extra_p: keep &= (u % extra_p != 0)
        cnt = cnt[keep]; u = u[keep]; cx = lut[u % q]
        rq = F[:, cx, :].transpose(1, 0, 2)[:, :, None, :] + F[cx][:, None, :, :] - F[None, :, :, :]
        n = cnt.sum((1, 2)); m1 = np.einsum('ulr,ulrd->ud', cnt, rq)/n[:, None]; m2 = np.einsum('ulr,ulrd->u', cnt, rq**2)/n
        noise = (m2 - (m1**2).sum(1))/(n-1)
        return float(((m1**2).sum(1) - noise).mean())/(C-1)**2
    base = sum(np.load('out2/v_%d.npz' % x)['c_%d_1' % q] for x in WIN); V0 = V_of(base, 210)
    # coarser: classes mod q only, mod 30, mod 210
    res = {'V210': V0, 'inc': {}}
    for p in PL:
        c = sum(np.load('out2/v_%d.npz' % x)['c_%d_%d' % (q, p)] for x in WIN)
        res['inc'][p] = V_of(c, 210*p, p) - V0
    cum = V0 + np.cumsum([res['inc'][p] for p in PL])
    out[q] = res
    print('q=%d V(210)=%.4f  inc(p): %s' % (q, V0, ' '.join('%d:%.4f' % (p, res['inc'][p]) for p in PL)))
    print('     additive estimate V(210*prod p<=P): ' + ' '.join('%d:%.3f' % (p, v) for p, v in zip(PL, cum)))
    print('     check: inc(11)+inc(13) = %.4f vs V(30030)-V(210) (direct) ; V(210)+inc11+inc13+inc17 = %.4f' % (res['inc'][11]+res['inc'][13], V0+res['inc'][11]+res['inc'][13]+res['inc'][17]))
json.dump(out, open('incr.json', 'w'), indent=1)
