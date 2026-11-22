"""First-order class chain mod 30030 (per-class gap law from review23/kth.py at the window's lambda), its last-digit
table for k = 1..3 conditioned on total distance <= DMAX, against the recounted data."""
import sys, math, numpy as np
gt, cf, DMAX = sys.argv[1], sys.argv[2], int(sys.argv[3])
PP = 30030
CLS = np.array([a for a in range(PP) if math.gcd(a, PP) == 1], np.int64); IDX = np.full(PP, -1); IDX[CLS] = np.arange(len(CLS))
G = np.load(gt); G = G / G.sum(1, keepdims=True)          # gap law per class, D = 2..DMAX, normalized
nD = G.shape[1]; Ds = 2 * np.arange(1, nD + 1)
dig = np.full(10, -1); dig[[1, 3, 7, 9]] = [0, 1, 2, 3]; cd = dig[CLS % 10]
S = np.zeros((4, len(CLS), nD + 1)); S[cd, np.arange(len(CLS)), 0] = 1.0   # state: start digit, class, distance/2
cnt = np.load(cf)['cnt']
for k in (1, 2, 3):
    new = np.zeros_like(S)
    for j, D in enumerate(Ds):
        h = D // 2; tgt = IDX[(CLS + D) % PP]; ok = (tgt >= 0) & (G[:, j] > 0)
        new[:, tgt[ok], h:] += S[:, ok, :nD + 1 - h] * G[ok, j][None, :, None]
    S = new
    tab = np.zeros((4, 4))
    for b in range(4): tab[:, b] = S[:, cd == b, :].sum((1, 2))
    law = tab / tab.sum(1, keepdims=True)
    tot = cnt[k-1].sum(0).reshape(4, 4).astype(float); N = tot.sum(1)[:, None]; dat = tot / N
    z = (law - dat) / np.sqrt(dat * (1 - dat) / N)
    print('k', k, 'chain max|z| %.1f, sum z^2 %.0f' % (np.abs(z).max(), (z**2).sum()), '1->1 chain %.4f data %.4f' % (law[0, 0], dat[0, 0]))
