"""diagnostic: pooled (noise-corrected at h=0) against cross-part A(h), q=5 k=1, response stretch 5e9"""
import sys, math, numpy as np
sys.path.insert(0, '.')
import map as M
S = M.MSet('new=val_new'); q, k = 5, 1
G = M.group(q); C = G.C; B, dg, P, pi = S.basis(q, k); D = B.shape[1]
F = M.yvec(P, pi) @ B; Fm = F.reshape(C*C, D)
Q = 510510; RU, _ = M.rough_index(Q, 17); NR = len(RU); cxu = G.lut[RU % q]
parts = 4
acc = {'n': np.zeros((parts, NR), np.int64), 1: np.zeros((parts, NR*C*C), np.int32)}
M.count_stretch(5*10**9, q, [1], parts, acc)
Sp = [M.signed_from_counts(acc[1][p].reshape(NR, C, C), cxu, C).reshape(NR, C*C) for p in range(parts)]
n = acc['n'].astype(float)
Stot = sum(Sp); ntot = n.sum(0)
for d in (0, 2):
    m1 = (Stot @ Fm[:, d]) / ntot
    full = np.zeros(Q); full[RU] = m1; Ff = np.fft.rfft(full); Apool = np.fft.irfft((Ff*Ff.conj()).real, n=Q)/NR
    Fp = []
    for p in range(parts):
        full = np.zeros(Q); full[RU] = (Sp[p] @ Fm[:, d]) / n[p]; Fp.append(np.fft.rfft(full))
    sF = sum(Fp); s2 = sum((f*f.conj()).real for f in Fp)
    Across = np.fft.irfft(((sF*sF.conj()).real - s2)/(parts*(parts-1)), n=Q)/NR
    diff = Apool - Across
    print('dir', d, 'A(0) pooled %.5f cross %.5f' % (Apool[0], Across[0]))
    for hs in (1, 2, 4, 6, 10, 20, 30, 60, 100, 200, 1000, 10000, 100000):
        print('   h=%6d  pooled %+.6f cross %+.6f diff %+.6f' % (hs, Apool[hs], Across[hs], diff[hs]))
    # where does the difference sit: sums of diff over ranges
    for lo, hi in ((1, 100), (100, 1000), (1000, 10000), (10000, 100000), (100000, Q//2)):
        print('   sum diff over h in [%d,%d): %+.5f ; mean %+.2e' % (lo, hi, diff[lo:hi].sum(), diff[lo:hi].mean()))
