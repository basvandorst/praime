"""Review 26: first-order chain on classes mod 30030 from an own per-class gap table T[class, D] (mylaw.py z^0 term).
Forward: E_k[c_end, D] = P(k-th neighbour at distance D, in class c_end), start class uniform. Paths within DMAX.
usage: chain.py table.npz lam_index KMAX out.npy"""
import sys, math
import numpy as np
Z = np.load(sys.argv[1]); li = int(sys.argv[2]); KMAX = int(sys.argv[3])
T = Z['T'][li].astype(np.float64); Ds = Z['D']; ND = len(Ds); PP = 30030
CLS = np.array([a for a in range(PP) if math.gcd(a, PP) == 1]); idx = np.full(PP, -1); idx[CLS] = np.arange(len(CLS))
nidx = idx[(CLS[:, None] + Ds[None, :]) % PP]
assert np.all(T[nidx < 0] == 0)
E = np.zeros((len(CLS), ND))
for j in range(ND):
    m = nidx[:, j] >= 0; np.add.at(E[:, j], nidx[m, j], T[m, j] / len(CLS))
out = np.zeros((KMAX, ND)); out[0] = E.sum(0)
for k in range(2, KMAX + 1):
    E2 = np.zeros_like(E)
    for j in range(ND - 1):
        if not E[:, j].any():
            continue
        for j2 in range(ND - j - 1):
            w = E[:, j] * T[:, j2]; m = (w > 0) & (nidx[:, j2] >= 0)
            np.add.at(E2[:, j + j2 + 1], nidx[m, j2], w[m])
    E = E2; out[k - 1] = E.sum(0)
np.save(sys.argv[4], out)
print('done', out.sum(1))
