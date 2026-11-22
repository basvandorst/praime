"""Two-step return P(D_2 = 0 mod q) of round 74's theory chain (exact sum over the table, uniform start classes)."""
import math, numpy as np
TH = np.load('r72_table.npy'); TH = TH / TH.sum(1, keepdims=True); NG = TH.shape[1]; gaps = np.arange(2, 2 * NG + 1, 2)
Q = 30030; CLS = np.array([a for a in range(Q) if math.gcd(a, Q) == 1]); IDX = np.full(Q, -1); IDX[CLS] = np.arange(len(CLS))
for q in (3, 5, 7, 11, 13):
    F = np.stack([TH[:, gaps % q == r].sum(1) for r in range(q)], 1)      # forward law of class b, by gap mod q
    tot = 0.0
    for i, c in enumerate(CLS):
        w = TH[IDX[(-c) % Q]]                                           # backward gap g1 from c
        nz = w > 0; g1 = gaps[nz]; c2 = (c - g1) % Q
        tot += float(np.sum(w[nz] * F[IDX[(-c2) % Q], (-g1) % q]))
    print(q, round(tot / len(CLS), 5))
