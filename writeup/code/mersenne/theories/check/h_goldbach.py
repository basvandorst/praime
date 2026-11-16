# H: exact means via counting pairs with p+q in the window (prefix counts), plus a random-sample direct count
import numpy as np
from math import isqrt

LO, HI = 2_000_000, 4_000_000
isp = np.ones(HI + 1, dtype=bool)
isp[:2] = False
for i in range(2, isqrt(HI) + 1):
    if isp[i]:
        isp[i * i::i] = False
pr = np.nonzero(isp)[0]
P = {a: pr[(pr % 12 == a) & (pr >= 5)] for a in (1, 5, 7, 11)}


def pairsum(a, b):
    # number of ordered (p, q), p = a, q = b mod 12, with LO <= p + q <= HI
    p = P[a]
    q = P[b]
    hi = np.searchsorted(q, HI - p, side="right")
    lo = np.searchsorted(q, LO - p, side="left")
    return int(np.sum(np.maximum(hi - lo, 0)))


ns = np.arange(LO, HI + 1)
ncount = {c: int(np.count_nonzero(ns % 12 == c)) for c in (0, 2, 6, 10)}
D = {
    0: pairsum(1, 11) + pairsum(11, 1) - pairsum(5, 7) - pairsum(7, 5),
    2: pairsum(1, 1) - pairsum(7, 7),
    6: pairsum(1, 5) + pairsum(5, 1) - pairsum(7, 11) - pairsum(11, 7),
    10: pairsum(5, 5) - pairsum(11, 11),
}
# S: sum over n of #{s prime >= 5, s^2 < n, n - s^2 prime}, split by n mod 12
Ssum = {c: 0 for c in (0, 2, 6, 10)}
for s in pr[(pr >= 5) & (pr * pr < HI)]:
    s2 = int(s) * int(s)
    lo = max(LO - s2, 2)
    hi = HI - s2
    if hi < lo:
        continue
    r = pr[(pr >= lo) & (pr <= hi)]
    n = r + s2
    for c in Ssum:
        Ssum[c] += int(np.count_nonzero(n % 12 == c))
print("number of n per class", ncount)
for c in (0, 2, 6, 10):
    print(f"class {c}: exact mean Delta = {D[c] / ncount[c]:.4f}, exact mean S = {Ssum[c] / ncount[c]:.4f}")

# random sample, direct per-n count
rng = np.random.default_rng(12345)
cls = {1: 1, 5: 5, 7: 7, 11: 11}
pm = pr[pr >= 5]
pm12 = pm % 12
sq = pm[pm * pm < HI]
for c in (0, 2, 6, 10):
    cand = ns[ns % 12 == c]
    samp = rng.choice(cand, 300, replace=False)
    dv = []
    sv = []
    for n in samp:
        n = int(n)
        ok = pm < n
        pp = pm[ok]
        m = isp[n - pp]
        cl = pm12[ok][m]
        cnt = {a: int(np.count_nonzero(cl == a)) for a in (1, 5, 7, 11)}
        # cnt[a] = ordered pairs with p = a (q = n - p class determined)
        if c == 0:
            d = cnt[1] + cnt[11] - cnt[5] - cnt[7]
        elif c == 2:
            d = cnt[1] - cnt[7]
        elif c == 6:
            d = cnt[1] + cnt[5] - cnt[7] - cnt[11]
        else:
            d = cnt[5] - cnt[11]
        # careful: q must be >= 5 too (q = 2 or 3 not in classes); q = n - p >= 5 since n large and q prime of class
        dv.append(d)
        ss = sq[sq * sq < n]
        sv.append(int(np.count_nonzero(isp[n - ss * ss])))
    dv = np.array(dv, float)
    sv = np.array(sv, float)
    print(f"sample class {c}: Delta {dv.mean():.2f} +- {dv.std(ddof=1) / np.sqrt(300):.2f}, "
          f"S {sv.mean():.2f} +- {sv.std(ddof=1) / np.sqrt(300):.2f}")
