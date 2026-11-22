"""Independent check of round 31: closed form of S4(d)=S({0,2,d,d+2}), the variance law, and the bias of
sliding-window variance estimates (demeaned within blocks) under the law itself."""
import math, json, sys
import numpy as np
from lib import primes_upto, sing

S2 = 1.3203236316937391  # 2 * C2
PR = primes_upto(20_000_000).astype(np.float64)
P5 = PR[PR >= 5]
G4 = 8 * 27 / 16 * math.exp(float(np.sum(np.log1p(-4 / P5) - 4 * np.log1p(-1 / P5))))
print('G4', G4)

DMAX = 19_000_010


def S4_closed(dmax):
    d = np.arange(dmax + 1, dtype=np.int64)
    S = np.where(d % 6 == 0, G4, 0.0); S[0] = 0.0
    del d
    for p in primes_upto(dmax + 3):
        p = int(p)
        if p < 5:
            continue
        S[0::p] *= (1 - 2 / p) / (1 - 4 / p)          # p | d
        f = (1 - 3 / p) / (1 - 4 / p)
        S[2::p] *= f                                   # p | d-2
        S[p - 2::p] *= f                               # p | d+2
    S[0] = 0.0
    return S


# 1) closed form vs direct series (own direct product to 2e6 + tail)
prs = primes_upto(2_000_000)
Ssmall = S4_closed(400)
worst = 0
for d in list(range(1, 61)) + [66, 90, 120, 150, 174, 210, 222, 330, 390]:
    direct = sing([0, 2, d, d + 2], prs=prs) if d != 2 else sing([0, 2, 4], prs=prs)
    c = Ssmall[d]
    if direct == 0 or c == 0:
        assert direct == 0 and c == 0, (d, direct, c)
    else:
        worst = max(worst, abs(c / direct - 1))
print('closed form vs direct: worst rel diff', worst)

mode = sys.argv[1] if len(sys.argv) > 1 else 'law'
S = S4_closed(DMAX)
A = np.cumsum(S)                                       # A[L-1] = sum_{d<L} S4(d)
dS = S * np.arange(len(S), dtype=np.float64)
del S
Bc = np.cumsum(dS); del dS


def cum(L):                                            # sum_{d<L} (L-d) S4(d)
    L = np.asarray(L)
    return L * A[L - 1] - Bc[L - 1]


def v(H, lam):
    r = S2 / lam ** 2
    return 1 - H * r + 2 * cum(H) / (H * r * lam ** 4)


def window_law(H, L, W, n=400):
    """average of v over [L, L+W) weighted by twin density 1/ln^2 x."""
    x = np.linspace(L, L + W, n + 1)
    lam = np.log(x); w = 1 / lam ** 2
    return float(np.sum(w * np.array([v(H, l) for l in lam])) / np.sum(w))


out = {'G4': G4, 'closed_vs_direct_worst': worst}
for name in ([] if mode == 'bias' else []): pass
for name, (L, W) in ({} if mode == 'bias' else {'train': (10 ** 7, 9 * 10 ** 7), 'A': (10 ** 10, 5 * 10 ** 7)}).items():
    for H in (1000, 10000, 100000, 1000000):
        th = [window_law(H, L + b * W // 10, W // 10, 40) for b in range(10)]
        out['%s_%d' % (name, H)] = float(np.mean(th))
        print(name, H, 'law %.4f' % np.mean(th), flush=True)
for H in (() if mode == 'bias' else (10 ** 6, 10 ** 7)):
    th = [window_law(H, 10 ** 9 + k * 2 * 10 ** 8, 2 * 10 ** 8, 40) for k in range(10)]
    out['1e9-3e9_%d' % H] = float(np.mean(th)); print('1e9-3e9', H, '%.4f' % np.mean(th))

# 2) bias of the sliding estimator used in r31.py (np.var of all sliding counts of one block, demeaned):
#    E[var] = V(H) - Var(mean of the N sliding counts), with Cov(T_H(0), T_H(s)) = (V(s+H) + V(|s-H|) - 2 V(s))/2
def Vfun(Ls, lam):
    r = S2 / lam ** 2
    Ls = np.asarray(Ls, dtype=np.int64)
    out = np.zeros(len(Ls))
    m = Ls > 0
    L = Ls[m].astype(np.float64)
    out[m] = L * r - L * L * r * r + 2 * cum(Ls[m]) / lam ** 4
    return out


for H, N in ((100000, 9 * 10 ** 6 - 100000 + 1), (1000000, 9 * 10 ** 6 + 1), (10 ** 7, 9 * 10 ** 6 + 1)):
    lam = math.log(5.5e7)                              # middle of train, one block
    VH = Vfun([H], lam)[0]
    tot = 0.0
    for a in range(1, N, 500000):
        s = np.arange(a, min(N, a + 500000), dtype=np.int64)
        C = 0.5 * (Vfun(s + H, lam) + Vfun(np.abs(s - H), lam) - 2 * Vfun(s, lam))
        tot += float(np.sum((N - s) * C))
    varmean = (N * VH + 2 * tot) / N ** 2
    mean = H * S2 / lam ** 2
    out['sliding_bias_%d' % H] = {'law_v': VH / mean, 'expected_sliding_v': (VH - varmean) / mean}
    print('sliding H=%d: law v %.4f, expected value of the sliding estimate %.4f' % (H, VH / mean, (VH - varmean) / mean), flush=True)
    del s, C
json.dump(out, open('/home/bas/scratch/study-review5/r31_%s.json' % mode, 'w'), indent=1)
