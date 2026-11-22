"""Round 31, addendum 2: twin-count variance in non-overlapping intervals of 1e6 and 1e7 in [1e9, 3e9)."""
import os, json, math
import numpy as np
from scipy.stats import t as student_t
import common as C
import theory2 as T2

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
L0, W0 = 10 ** 9, 2 * 10 ** 9
HS = [10 ** 6, 10 ** 7]
DMAX = max(HS) + 4
PR = C.small_primes(20_000_000).astype(np.float64)
G4 = 8 * 27 / 16 * float(np.exp(np.sum(np.log1p(-4 / PR[PR >= 5]) - 4 * np.log1p(-1 / PR[PR >= 5]))))
S4 = np.where(np.arange(DMAX + 1) % 6 == 0, G4, 0.0); S4[0] = 0.0
for p in C.small_primes(DMAX + 3):          # the closed form of r31.py (checked there against S_many for d <= 200)
    p = int(p)
    if p < 5:
        continue
    S4[p::p] *= (1 - 2 / p) / (1 - 4 / p)
    f = (1 - 3 / p) / (1 - 4 / p)
    S4[2::p] *= f
    S4[p - 2::p] *= f
S2 = T2.S1([0, 2])
cum = {H: float(np.sum((H - np.arange(1, H)) * S4[1:H])) for H in HS}


def law(H, L, W):
    return T2.window_avg(L, W, lambda lam: 1 - H * S2 / lam ** 2 + 2 * cum[H] / (H * (S2 / lam ** 2) * lam ** 4))
counts = {H: [] for H in HS}; mids = {H: [] for H in HS}
carry = None
for s in range(200):
    lo = L0 + s * 10 ** 7
    M = C.primes_mask(lo, 10 ** 7 + 2)
    tw = M[:-2] & M[2:]
    tw = tw[:10 ** 7]
    c6 = tw.reshape(10, 10 ** 6).sum(1)
    counts[10 ** 6] += c6.tolist(); mids[10 ** 6] += (lo + np.arange(10) * 10 ** 6 + 5e5).tolist()
    counts[10 ** 7].append(int(c6.sum())); mids[10 ** 7].append(lo + 5e6)
out = {}
for H in HS:
    c = np.array(counts[H], float); m = np.array(mids[H])
    e = H * S2 / np.log(m) ** 2
    res = c - e
    K = 10; per = len(c) // K
    meas = [float(np.var(res[k * per:(k + 1) * per], ddof=1) / e[k * per:(k + 1) * per].mean()) for k in range(K)]
    th = [law(H, L0 + k * W0 // K, W0 // K) for k in range(K)]
    r = np.array(meas) - np.array(th)
    t9 = float(r.mean() / (r.std(ddof=1) / math.sqrt(K)))
    out[str(H)] = {'intervals': len(c), 'meas_pooled': float(np.var(res, ddof=1) / e.mean()), 'meas': float(np.mean(meas)),
                   'law': float(np.mean(th)), 't9': t9, 'p': float(2 * student_t.sf(abs(t9), K - 1))}
    print(H, out[str(H)], flush=True)
json.dump(out, open(os.path.join(OUT, 'r31b.json'), 'w'), indent=1)
