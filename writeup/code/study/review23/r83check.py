"""Round 83 check with own data histogram (c81.json), own law (t81.npz), own Gamma and independent-gap comparators, and
an extra comparator: the first-order chain on classes mod 30030 (own per-class gap law, t81_gaptable.npy)."""
import json, math
import numpy as np
from scipy import stats

DMAX = 240; PP = 30030
c = json.load(open('c81.json')); t = np.load('t81.npz'); T = np.load('t81_gaptable.npy')
Ds = t['D']; lam = float(t['lam'])
CLS = np.array([a for a in range(PP) if math.gcd(a, PP) == 1]); IDX = np.full(PP, -1); IDX[CLS] = np.arange(len(CLS))
bands = [(lo, lo + 10) for lo in range(2, DMAX + 1, 12)]


def band(v):          # v indexed by even D = 2..DMAX
    return np.array([v[(Ds >= lo) & (Ds <= hi)].sum() for lo, hi in bands])


def chi2(cnt, e):     # pool forward until expected >= 50, remainder into the last kept band
    rows = []; a = b = 0.0
    for x, y in zip(cnt, e):
        a += x; b += y
        if b >= 50:
            rows.append([a, b]); a = b = 0.0
    if b > 0:
        rows[-1][0] += a; rows[-1][1] += b
    return sum((x - y) ** 2 / y for x, y in rows), len(rows) - 1


# class chain: V[class, D/2]
Tn = T / T.sum(1, keepdims=True)
V = np.zeros((len(CLS), DMAX // 2 + 1)); V[:, 0] = 1 / len(CLS)
chain = {}
for k in range(1, 5):
    W = np.zeros_like(V)
    for j in range(DMAX // 2):
        g = 2 * (j + 1); dst = IDX[(CLS + g) % PP]; ok = dst >= 0
        contrib = V[ok, :DMAX // 2 + 1 - (j + 1)] * Tn[ok, j][:, None]
        np.add.at(W, (dst[ok][:, None], np.arange(j + 1, DMAX // 2 + 1)[None, :]), contrib)
    V = W; chain[k] = V.sum(0)[1:]
one = t['dist'][0] / t['dist'][0].sum()
for k in range(1, 5):
    h = np.array(c['res'][str(k)]['hist_even_D'])[1:]      # D = 2..DMAX
    N = h.sum(); cnt = band(h.astype(float))
    out = {}
    law = t['dist'][k - 1]; out['law'] = chi2(cnt, N * band(law / law.sum()))
    gk = stats.gamma.pdf(Ds, k, scale=lam); out['gamma'] = chi2(cnt, N * band(gk / gk.sum()))
    # independent gaps by explicit enumeration of the k-fold sum (own indexing)
    f = np.zeros(DMAX // 2 + 1); f[0] = 1
    for _ in range(k):
        f2 = np.zeros_like(f)
        for j in range(DMAX // 2):
            f2[j + 1:] += f[:DMAX // 2 - j] * one[j]
        f = f2
    ind = f[1:]; out['indep'] = chi2(cnt, N * band(ind / ind.sum()))
    ch = chain[k]; out['class_chain'] = chi2(cnt, N * band(ch / ch.sum()))
    # repeat probabilities of the chain, for the record
    rep = {q: float(ch[Ds % q == 0].sum() / ch.sum()) for q in (3, 5, 7, 11, 13)}
    print(k, 'N', N, {n: (round(x, 1), d) for n, (x, d) in out.items()}, 'chain rep', {q: round(v, 4) for q, v in rep.items()})
