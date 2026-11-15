#!/usr/bin/env python3
"""Transfer: the six frozen classes on the Mersenne exponents.

Universe: primes 67 <= k < 83,880,077 (2^k - 1 >= 2^64, all exponents verified by GIMPS).
Null (Wagstaff, as in MERSENNE.md section 3): w = e^gamma log2(a k)/k, a = 2 (k = 3 mod 4), 6 (k = 1 mod 4).
E for H1, H3, H4, H5: exact over all primes. H2, H6: exact below 10^6, above by a sample of
primes drawn with probability proportional to w (importance sample).
Also the 52nd exponent 136,279,841 (beyond the verified bound): in class or not.
"""
import math, os, json
import numpy as np
import features as FT
from scipy.stats import poisson

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253, 4423,
       9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243, 110503, 132049, 216091, 756839, 859433,
       1257787, 1398269, 2976221, 3021377, 6972593, 13466917, 20996011, 24036583, 25964951,
       30402457, 32582657, 37156667, 42643801, 43112609, 57885161, 74207281, 77232917, 82589933]
NEW = 136279841
LIM = 83880077
EG = math.exp(0.5772156649015329)


def primes_upto(n):
    s = np.ones(n // 2, bool); s[0] = False          # s[i] <-> 2i+1
    for i in range(3, int(n ** 0.5) + 1, 2):
        if s[i // 2]:
            s[i * i // 2::i] = False
    return 2 * np.flatnonzero(s) + 1


def glide_vec(k):
    v = k.copy(); g = np.zeros(len(k), np.int64); act = np.arange(len(k)); t = 0
    while len(act):
        t += 1
        x = v[act]
        x = np.where(x & 1, 3 * x + 1, x >> 1)
        v[act] = x
        d = x < k[act]
        g[act[d]] = t
        act = act[~d]
    return g


def cf_vec(k):
    a0 = np.floor(np.sqrt(k.astype(np.float64))).astype(np.int64)
    a0 = np.where((a0 + 1) ** 2 <= k, a0 + 1, a0); a0 = np.where(a0 ** 2 > k, a0 - 1, a0)
    m = np.zeros(len(k), np.int64); d = np.ones(len(k), np.int64); a = a0.copy(); mx = np.zeros(len(k), np.int64)
    act = np.arange(len(k))
    while len(act):
        mm, dd, aa, kk, A0 = m[act], d[act], a[act], k[act], a0[act]
        mm = dd * aa - mm; dd = (kk - mm * mm) // dd; aa = (A0 + mm) // dd
        m[act], d[act], a[act] = mm, dd, aa
        fin = aa == 2 * A0
        mx[act[~fin]] = np.maximum(mx[act[~fin]], aa[~fin])
        act = act[~fin]
    return mx / np.sqrt(k.astype(np.float64))


def main():
    P = primes_upto(LIM)
    P = P[P >= 67].astype(np.int64)
    a = np.where(P % 4 == 3, 2.0, 6.0)
    w = EG * np.log2(a * P) / P
    ex = np.array([e for e in EXP if 67 <= e < LIM], np.int64)
    isexp = np.isin(P, ex)
    print(f"universe: {len(P)} primes in [67, {LIM}), {isexp.sum()} exponents, E_total = {w.sum():.2f}")
    rng = np.random.default_rng(5)
    lo = P < 10 ** 6
    hiidx = np.flatnonzero(~lo)
    samp = rng.choice(hiidx, size=60000, p=w[hiidx] / w[hiidx].sum())
    Whi = w[hiidx].sum()
    res = {}
    for name, fn, test in FT.FROZEN:
        key = name.split()[0]
        if key in ("H1", "H3", "H4", "H5"):
            if key == "H5":
                v = glide_vec(P).astype(float)
            else:
                v = fn(P)
            cl = test(v)
            E = w[cl].sum(); O = int((cl & isexp).sum())
        else:
            if key == "H6":
                vlo = cf_vec(P[lo]); vs = cf_vec(P[samp]); vex = cf_vec(ex)
            else:
                vlo = fn(P[lo]); vs = fn(P[samp]); vex = fn(ex)
            cllo = np.where(np.isnan(vlo), False, test(np.nan_to_num(vlo, nan=-9)))
            cls = np.where(np.isnan(vs), False, test(np.nan_to_num(vs, nan=-9)))
            clex = np.where(np.isnan(vex), False, test(np.nan_to_num(vex, nan=-9)))
            E = w[lo][cllo].sum() + Whi * cls.mean()
            O = int(clex.sum())
        vnew = fn(np.array([NEW], np.int64))[0] if key != "H5" else float(glide_vec(np.array([NEW], np.int64))[0])
        innew = bool(test(np.array([vnew]))[0]) if not np.isnan(vnew) else False
        p_up = poisson.sf(O - 1, E)
        res[name] = dict(O=O, E=float(E), p_upper=float(p_up), new_value=float(vnew), new_in_class=innew)
        print(f"{name:42s} O {O:2d}  E {E:6.2f}  O/E {O / E:5.2f}  p_up {p_up:.3g}   52nd: value {vnew:.4f} in class {innew}")
    json.dump(res, open(os.path.join(HERE, "mersenne.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
