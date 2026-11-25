# Independent constants for level-1 forms: c = B(1)/L(1,Sym^2 f), L(1,Sym^2 f) = 2^{2k-1} pi^{k+1} <f,f> / Gamma(k)
# (Rankin; <f,f> = int_{SL2(Z)\H} |f|^2 y^k dx dy/y^2), B_p = (1-e X)/((1-X)(1-(lam^2-2)X+X^2)), e = lam^2 - 1, X = 1/p.
# Second route: direct conditionally convergent product prod_{p<=P} (1 - (lam_p^2-1)/p) with log-Cesaro mean.
import os, sys, time, math, numpy as np, cypari2, json
pari = cypari2.Pari(); pari.allocatemem(800 * 10 ** 6)
pari.set_real_precision(30)
P = 200000
def primes(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]
pr = primes(P)
res = {}
for k in (12, 16, 18, 20, 22, 26):
    t0 = time.time()
    mf = pari('mfinit([1,%d],1)' % k)
    f = pari.mfbasis(mf)[0]
    co = pari.mfcoefs(f, P)
    lam = np.array([float(co[int(p)]) / float(p) ** ((k - 1) / 2) for p in pr])
    X = 1.0 / pr
    e = lam ** 2 - 1
    Bp = (1 - e * X) / ((1 - X) * (1 - (lam ** 2 - 2) * X + X ** 2))
    logB = np.cumsum(np.log(Bp))
    B = math.exp(logB[-1])
    pet = float(pari.mfpetersson(pari.mfsymbol(mf, f)))
    LS = 2 ** (2 * k - 1) * math.pi ** (k + 1) * pet / math.gamma(k)
    c = B / LS
    # direct product
    lf = np.log(1 - e * X)
    cum = np.cumsum(lf)
    # log-Cesaro mean over P in [P/10, P] with weights dP/P, then bias correction exp(+nu/(sqrt P ln P)), nu = 1
    m = pr >= P / 10
    w = 1.0 / pr[m]
    direct_end = math.exp(cum[-1])
    direct_smooth = math.exp((cum[m] * w).sum() / w.sum())
    nu = 1.0
    bias = np.exp(nu / (np.sqrt(pr[m]) * np.log(pr[m])))
    direct_smooth_b = math.exp(((cum[m] + np.log(bias)) * w).sum() / w.sum())
    res[k] = dict(B=B, pet=pet, LS1=LS, c=c, direct_end=direct_end, direct_smooth=direct_smooth, direct_smooth_bias=direct_smooth_b, lam2=lam[0])
    print('k=%d lam_2=%+.4f <f,f>=%.10e L(1,Sym2)=%.10f B(1)=%.8f c=%.6f | direct P=%d: end %.5f, log-Cesaro %.5f, +bias %.5f  (%.0fs)'
          % (k, lam[0], pet, LS, B, c, P, direct_end, direct_smooth, direct_smooth_b, time.time() - t0))
json.dump(res, open(os.path.dirname(os.path.abspath(__file__))+'/r6a_level1.json', 'w'), indent=1)
