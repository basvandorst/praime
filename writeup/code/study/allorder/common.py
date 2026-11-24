"""Shared helpers. The within basis (character pairs, Loewdin-orthonormalised inside the exact within projector) is the
construction of review44/dirs.py, re-typed here so that per-direction numbers are comparable."""
import numpy as np, math, cmath
STR = {'A': 150_000_000_000, 'B': 5_000_000_000}
SPAN = 100_000_000
GEN = {3: 2, 5: 2, 7: 3, 9: 2, 11: 2, 13: 2}
def red(q): return [a for a in range(q) if math.gcd(a, q) == 1]
def chars(q):
    r = red(q); C = len(r); g = GEN[q]; lg = {}; x = 1
    for j in range(C): lg[x] = j; x = x * g % q
    return np.array([[cmath.exp(2j * math.pi * k * lg[a] / C) for a in r] for k in range(C)])
def projector(P, pi):
    C = len(pi); R = np.zeros((C * C, C * C)); V = np.zeros((C * C, C))
    for a in range(C):
        s = np.sqrt(P[a]); R[a*C:(a+1)*C, a*C:(a+1)*C] = np.eye(C) - np.outer(s, s)
        for b in range(C): V[a*C+b, b] = math.sqrt(pi[a] * P[a, b])
    return R - R @ V @ np.linalg.pinv(V.T @ R @ V) @ V.T @ R
def basis(q, P, pi):
    """columns: orthonormal within directions; labels '(k,l)type'; parity under (a,b)->(-a,-b) = chi_k(-1) chi_l(-1)"""
    C = len(pi); ch = chars(q); PW = projector(P, pi); w = np.sqrt(pi[:, None] * P).ravel()
    r = red(q); im1 = r.index(q - 1)
    vecs, labels, par = [], [], []; done = set()
    for k in range(1, C):
        for l in range(1, C):
            if (k, l) in done: continue
            kc, lc = (C - k) % C, (C - l) % C; done.add((k, l)); done.add((kc, lc))
            F = np.outer(ch[k], np.conj(ch[l]))
            typ = 'ratio' if k == l else ('product' if (k + l) % C == 0 else 'mixed')
            lab = '(%d,%d)%s' % (k, l, typ)
            sg = int(round((ch[k][im1] * np.conj(ch[l][im1])).real))
            parts = [F.real] if (k, l) == (kc, lc) else [F.real, F.imag]
            for Fp in parts: vecs.append(PW @ (Fp.ravel() * w)); labels.append(lab); par.append(sg)
    Vm = np.array(vecs).T
    G = Vm.T @ Vm; e, U = np.linalg.eigh(G); S = U @ np.diag(e ** -0.5) @ U.T
    B = Vm @ S
    assert np.allclose(B.T @ B, np.eye(B.shape[1]), atol=1e-8) and B.shape[1] == (C - 1) ** 2
    assert np.allclose(B @ B.T, PW, atol=1e-6)
    return B, labels, np.array(par)
def yvec(P, pi):
    """Y[a, b] = within-space increment vector of one step a -> b"""
    C = len(pi); Y = np.zeros((C, C, C * C))
    for a in range(C):
        for b in range(C):
            for b2 in range(C): Y[a, b, a*C+b2] = ((1.0 if b2 == b else 0.0) - P[a, b2]) / math.sqrt(pi[a] * P[a, b2])
    return Y
def load(name):
    return np.load('primes_%s.npy' % name)
