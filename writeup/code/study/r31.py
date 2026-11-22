"""Round 31 (preregister/round31.txt): the Montgomery-Soundararajan analogue for twin primes."""
import os, json, math
import numpy as np
from scipy.stats import t as student_t
import theory2 as T2
import common as C

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
HS = [1000, 10000, 100000]
DMAX = max(HS) + 4
PR = C.small_primes(20_000_000).astype(np.float64)
G4 = 8 * 27 / 16 * float(np.exp(np.sum(np.log1p(-4 / PR[PR >= 5]) - 4 * np.log1p(-1 / PR[PR >= 5]))))


def S4_closed():
    d = np.arange(DMAX + 1)
    S = np.where(d % 6 == 0, G4, 0.0)
    S[0] = 0.0
    for p in C.small_primes(DMAX + 3):
        p = int(p)
        if p < 5:
            continue
        S[p::p] *= (1 - 2 / p) / (1 - 4 / p)                     # p | d
        f = (1 - 3 / p) / (1 - 4 / p)
        S[2 % p::p] *= f if p != 2 else 1                         # p | d - 2
        S[(p - 2) % p::p] *= f                                    # p | d + 2
    return S


S4 = S4_closed()
# check against the direct singular series
for d in range(6, 201, 6):
    direct = T2.S1([0, 2, d, d + 2])
    assert abs(direct - S4[d]) < 1e-6 * direct, (d, direct, S4[d])
S2 = T2.S1([0, 2])
cum = {}
for H in HS:
    dd = np.arange(1, H)
    cum[H] = float(np.sum((H - dd) * S4[1:H]))


def law(H, L, W):
    def f(lam):
        r = S2 / lam ** 2
        return 1 - H * r + 2 * cum[H] / (H * r * lam ** 4)
    return T2.window_avg(L, W, f)


out = {'G4': G4}
for win in ('train', 'A'):
    L, W = C.WINDOWS[win]
    M = C.primes_mask(L, W + 2)
    tw = (M[:-2] & M[2:]).astype(np.int8)                         # twin starting at position i
    res = {}
    for H in HS:
        meas, th = [], []
        for b in range(10):
            lo, n = b * W // 10, W // 10
            c = np.concatenate((np.zeros(1, np.int32), np.cumsum(tw[lo:lo + n], dtype=np.int32)))
            cnt = (c[H:] - c[:-H]).astype(np.float64)
            # detrend by the law's mean at each position (slow drift of lam inside the block)
            pos = L + lo + np.arange(len(cnt)) + H / 2
            e = H * S2 / np.log(pos) ** 2
            meas.append(float(np.var(cnt - e) / e.mean())); th.append(law(H, L + lo, n))
        r = np.array(meas) - np.array(th)
        t9 = float(r.mean() / (r.std(ddof=1) / math.sqrt(10)))
        res[str(H)] = {'meas': float(np.mean(meas)), 'law': float(np.mean(th)), 't9': t9, 'p': float(2 * student_t.sf(abs(t9), 9))}
        print(win, H, res[str(H)], flush=True)
    out[win] = res
json.dump(out, open(os.path.join(OUT, 'r31.json'), 'w'), indent=1)
