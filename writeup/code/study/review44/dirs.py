"""Direction-resolved within part: per character-pair group, the measured block mean square at H = 1e4, 1e5, 1e6 and its
slope, against the linear-response slope -V_g/lam per direction (V_g from the empirical responses mod 510510)."""
import numpy as np, math, json, cmath
WIN = [6*10**9 + i*10**8 for i in range(24)]
lam = np.mean([math.log(x+5e7) for x in WIN])
QR = 510510
okR = np.ones(QR, bool)
for p in (2,3,5,7,11,13,17): okR[::p] = False
RU = np.flatnonzero(okR); NR = len(RU)
def red(q): return [a for a in range(q) if math.gcd(a, q) == 1]
GEN = {3: 2, 5: 2, 7: 3, 10: 3}
def chars(q):
    r = red(q); C = len(r); g = GEN[q]; lg = {}
    x = 1
    for j in range(C): lg[x] = j; x = x*g % q
    return np.array([[cmath.exp(2j*math.pi*k*lg[a]/C) for a in r] for k in range(C)])   # chi_k(a)
def projectors(P, pi):
    C = len(pi); R = np.zeros((C*C, C*C)); V = np.zeros((C*C, C))
    for a in range(C):
        s = np.sqrt(P[a]); R[a*C:(a+1)*C, a*C:(a+1)*C] = np.eye(C) - np.outer(s, s)
        for b in range(C): V[a*C+b, b] = math.sqrt(pi[a]*P[a, b])
    return R - R @ V @ np.linalg.pinv(V.T @ R @ V) @ V.T @ R
def basis(q, P, pi):
    """orthonormal basis of the within space from character pairs (Loewdin), with group labels"""
    C = len(pi); ch = chars(q); PW = projectors(P, pi); w = np.sqrt(pi[:, None]*P).ravel()
    vecs, labels = [], []
    done = set()
    for k in range(1, C):
        for l in range(1, C):
            if (k, l) in done: continue
            kc, lc = (C-k) % C, (C-l) % C; done.add((k, l)); done.add((kc, lc))
            F = np.outer(ch[k], np.conj(ch[l]))
            typ = 'ratio' if k == l else ('product' if (k+l) % C == 0 else 'mixed')
            lab = '(%d,%d)%s' % (k, l, typ)
            parts = [F.real] if (k, l) == (kc, lc) else [F.real, F.imag]
            for Fp in parts: vecs.append(PW @ (Fp.ravel()*w)); labels.append(lab)
    Vm = np.array(vecs).T                       # (C*C, D)
    G = Vm.T @ Vm; e, U = np.linalg.eigh(G); S = U @ np.diag(e**-0.5) @ U.T
    B = Vm @ S                                   # Loewdin-orthonormal, same span
    assert np.allclose(B.T @ B, np.eye(B.shape[1]), atol=1e-8) and B.shape[1] == (C-1)**2
    assert np.allclose(B @ B.T, PW, atol=1e-6)
    return B, labels
def yvec(P, pi):
    C = len(pi); Y = np.zeros((C, C, C*C))
    for a in range(C):
        for b in range(C):
            for b2 in range(C): Y[a, b, a*C+b2] = ((1.0 if b2 == b else 0.0) - P[a, b2])/math.sqrt(pi[a]*P[a, b2])
    return Y
out = {}
for q, k in ((5, 1), (7, 1), (3, 1), (10, 2)):
    r = red(q); C = len(r); key = lambda H: 'X_%d_%d_%d' % (q, k, H)
    Ttot = sum(np.load('out3/t_%d.npz' % x)[key(10**4)].sum(0) for x in WIN).astype(float)
    Pg = Ttot/Ttot.sum(1, keepdims=True); pig = Ttot.sum(1)/Ttot.sum()
    B, labels = basis(q, Pg, pig); groups = sorted(set(labels), key=labels.index)
    meas = {}
    for H in (10**4, 10**5, 10**6):
        acc = []   # per window, per direction mean square
        for x in WIN:
            X = np.load('out3/t_%d.npz' % x)[key(H)].astype(float); NB = X.shape[0]
            T = X.sum(0); P = T/T.sum(1, keepdims=True)
            Na = X.sum(2); E = Na[:, :, None]*P[None]
            y = ((X - E)/np.sqrt(E)).reshape(NB, C*C)
            z = y @ B                                     # (NB, D)
            acc.append((z**2).mean(0)*NB/(NB-1))
        acc = np.array(acc)                               # (24, D)
        meas[H] = acc
    res = {}
    # linear-response V per direction (k = 1 only, response data mod 510510)
    Vd = None
    if k == 1:
        lut = -np.ones(q, np.int64); lut[r] = np.arange(C)
        F = yvec(Pg, pig) @ B
        cnt = sum(np.load('out/w_%d.npz' % x)['rc_comp_%d' % q] + np.load('out/w_%d.npz' % x)['rc_prime_%d' % q] for x in WIN).reshape(NR, C, C).astype(float)
        cx = lut[RU % q]
        rq = F[:, cx, :].transpose(1, 0, 2)[:, :, None, :] + F[cx][:, None, :, :] - F[None, :, :, :]
        n = cnt.sum((1, 2)); m1 = np.einsum('ulr,ulrd->ud', cnt, rq)/n[:, None]
        m2d = np.einsum('ulr,ulrd->ud', cnt, rq**2)/n[:, None]
        Vd = ((m1**2) - (m2d - m1**2)/(n[:, None]-1)).mean(0)
    print('q=%d k=%d  (%d within directions; groups of character pairs, gen %d)' % (q, k, (C-1)**2, GEN[q]))
    tot = {H: meas[H].sum(1).mean()/(C-1)**2 for H in meas}
    print('   total within/dof: %s' % ' '.join('%.4f' % tot[H] for H in (10**4, 10**5, 10**6)))
    for g in groups:
        idx = [i for i, l in enumerate(labels) if l == g]; nd = len(idx)
        v = {H: meas[H][:, idx].sum(1)/nd for H in meas}         # per window, per-direction average in the group
        m = {H: v[H].mean() for H in v}; e = {H: v[H].std(ddof=1)/math.sqrt(24) for H in v}
        d = (v[10**6] - v[10**4])/math.log(100); s, se = d.mean(), d.std(ddof=1)/math.sqrt(24)
        line = '   %-14s n=%d  W(1e4) %.4f+-%.4f  W(1e5) %.4f+-%.4f  W(1e6) %.4f+-%.4f  slope/unit lnH %+.4f+-%.4f' % (g, nd, m[10**4], e[10**4], m[10**5], e[10**5], m[10**6], e[10**6], s, se)
        if Vd is not None:
            vg = Vd[idx].mean(); line += '  | lin-resp (y=17, asympt.) -V/lam %+.4f  excess %+.4f (%.1f sig)' % (-vg/lam, s + vg/lam, (s + vg/lam)/se)
        print(line)
        res[g] = dict(n=nd, W={str(H): [m[H], e[H]] for H in m}, slope=[s, se], V=None if Vd is None else float(Vd[idx].mean()))
    out['%d_%d' % (q, k)] = res
json.dump(out, open('dirs.json', 'w'), indent=1)
