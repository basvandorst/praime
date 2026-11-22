"""finite-lam hidden-residue model of gap residues mod 30, own code. M = modulus/2 classes of even gaps (15: q=3,5;
105: q=3,5,7). Joint of k+1 consecutive gap classes ∝ prod w(c_i) * prod_q (q - nu_q), projected to classes mod 30."""
import math, sys, json
import numpy as np
def Hb(p):
    p = np.asarray(p, float).ravel(); p = p[p > 0]; p = p / p.sum(); return float(-(p * np.log2(p)).sum())
def mi(J):
    return Hb(J.sum(axis=-1)) + Hb(J.reshape(-1, 15).sum(axis=0)) - Hb(J)
def model(lam, k, qs):
    mod = 2 * int(np.prod(qs)); M = mod // 2
    g = np.arange(2, int(60 * lam) + mod, 2)
    w = np.bincount((g // 2) % M, weights=np.exp(-g / lam), minlength=M)
    n = k + 1
    J = np.zeros(15 ** n)
    rest = np.array(np.meshgrid(*([np.arange(M)] * (n - 1)), indexing='ij')).reshape(n - 1, -1)  # classes 2..n
    for c1 in range(M):
        C = np.vstack([np.full(rest.shape[1], c1), rest])           # n x L
        pos = np.vstack([np.zeros(C.shape[1], np.int64), np.cumsum(2 * C, axis=0)])   # n+1 positions
        f = np.prod(w[C], axis=0)
        for q in qs:
            R = np.sort(pos % q, axis=0)
            nu = 1 + np.count_nonzero(np.diff(R, axis=0), axis=0)
            f = f * (q - nu)
        code = np.zeros(C.shape[1], np.int64)
        for i in range(n):
            code = code * 15 + C[i] % 15
        J += np.bincount(code, weights=f, minlength=15 ** n)
    return mi(J.reshape((15,) * n))
out = {}
for lam in (float(np.mean(np.log(np.linspace(1e10, 1.1e10, 1001)))), 23.412105568345556):
    for k in (1, 2, 3, 4):
        a = model(lam, k, (3, 5)); out['%.4f_k%d_q35' % (lam, k)] = a
        print('lam %.4f k %d q3,5 %.5f' % (lam, k, a), flush=True)
    for k in (1, 2, 3):
        b = model(lam, k, (3, 5, 7)); out['%.4f_k%d_q357' % (lam, k)] = b
        print('lam %.4f k %d q3,5,7 %.5f' % (lam, k, b), flush=True)
json.dump(out, open('r62_model.json', 'w'), indent=1)
