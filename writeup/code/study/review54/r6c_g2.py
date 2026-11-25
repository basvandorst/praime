# Genus-2 constants: exact local factors p <= P0 (PARI hyperellcharpoly / lfuneuler), first order P0 < p <= P (own point count).
import os, sys, math, json, time, numpy as np, cypari2
pari = cypari2.Pari(); pari.allocatemem(300 * 10 ** 6)
OUT = os.path.dirname(os.path.abspath(__file__))+'/'
P0 = 1000; P = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
def primes(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]
pr = primes(P)
curves = {'G249': ([0, 1, 1], [1, 0, 0, 1], 2), 'G277': ([0, -1, 2, -2, 0, 1], [1], 1)}   # f coeffs (low->high), h coeffs, points at infinity rule
res = {}
for key, (fc, hc, _) in curves.items():
    t0 = time.time()
    fs = '+'.join('(%d)*x^%d' % (c, i) for i, c in enumerate(fc)); hs = '+'.join('(%d)*x^%d' % (c, i) for i, c in enumerate(hc))
    Lg = pari('lfungenus2([%s,%s])' % (fs, hs)); N = int(pari('lfunparams(lfungenus2([%s,%s]))[1]' % (fs, hs)))
    lt = np.zeros(len(pr)); Bt = np.ones(len(pr)); lam1 = np.zeros(len(pr))
    for i, p in enumerate(pr):
        p = int(p); x = 1.0 / p
        if p <= P0 or N % p == 0:
            if N % p == 0:
                eu = pari.lfuneuler(Lg, p); poly = pari.denominator(eu); num = pari.numerator(eu)
                cf = [float(pari.polcoef(poly, j)) / float(num) for j in range(int(pari.poldegree(poly)) + 1)]
            else:
                cp = pari('hyperellcharpoly([Mod(1,%d)*(%s), Mod(1,%d)*(%s)])' % (p, fs, p, hs))
                a1 = -int(pari.polcoef(cp, 3)); a2 = int(pari.polcoef(cp, 2))
                cf = [1, -a1, a2, -p * a1, p * p]
            b = np.array([cf[j] / p ** (j / 2) for j in range(len(cf))], complex)
            n = len(b) - 1
            lam = np.zeros(n + 1, complex); lam[0] = 1
            for k in range(1, n + 1): lam[k] = -sum(b[j] * lam[k - j] for j in range(1, k + 1))
            Fv = (np.conj(lam) * b * x ** np.arange(n + 1)).sum().real
            alpha = 1 / np.roots(b[::-1]) if n > 0 else []
            rs = 1.0 + 0j
            for a_ in alpha:
                for b_ in alpha: rs *= 1 / (1 - a_ * np.conj(b_) * x)
            lt[i] = math.log(Fv / (1 - x)); Bt[i] = (Fv * rs).real; lam1[i] = -b[1].real
            if N % p == 0: print('   %s bad p=%d analytic local poly %s, F_p(1/p)/(1-1/p) = %.6f, B_p = %.6f' % (key, p, np.round(b.real, 5), Fv / (1 - x), Bt[i]))
            if p in (101, 1009 if P0 > 1009 else 101):
                pass
        else:
            # own point count: y^2 + h y = f  <=>  (2y + h)^2 = h^2 + 4 f  (p odd)
            xs = np.arange(p, dtype=np.int64)
            def ev(c):
                v = np.zeros(p, np.int64)
                for cc in reversed(c): v = (v * xs + cc) % p
                return v
            hv = ev(hc); fv = ev(fc)
            disc = (hv * hv + 4 * fv) % p
            sq = np.zeros(p, np.int8); sq[(xs * xs) % p] = 1
            leg = np.where(disc == 0, 0, 2 * sq[disc].astype(np.int64) - 1)
            # points at infinity: deg f = 6 or deg h = 3 -> roots of Y^2 + h3 Y - f6 ; else 1
            h3 = hc[3] if len(hc) > 3 else 0; f6 = fc[6] if len(fc) > 6 else 0
            if h3 % p or f6 % p:
                dd = (h3 * h3 + 4 * f6) % p
                inf = 1 + (0 if dd == 0 else (1 if sq[dd] else -1))
            else:
                inf = 1
            Np = p + int(leg.sum()) + inf
            a1 = p + 1 - Np
            l1 = a1 / math.sqrt(p)
            lt[i] = math.log(1 - (l1 * l1 - 1) * x); Bt[i] = 1.0; lam1[i] = l1
    # spot check of the own count against PARI at two primes
    for p in (1009, 1013):
        cp = pari('hyperellcharpoly([Mod(1,%d)*(%s), Mod(1,%d)*(%s)])' % (p, fs, p, hs))
        j = int(np.nonzero(pr == p)[0][0])
        print('   check p=%d: PARI a1 = %d, table lambda*sqrt p = %.3f' % (p, -int(pari.polcoef(cp, 3)), lam1[j] * math.sqrt(p)))
    cum = np.cumsum(lt); m = pr >= P / 10; w = 1.0 / pr[m]
    d_end = math.exp(cum[-1]); d_sm = math.exp((cum[m] * w).sum() / w.sum())
    d_smb = math.exp(((cum[m] + 2 / (np.sqrt(pr[m]) * np.log(pr[m]))) * w).sum() / w.sum())
    B = math.exp(np.log(Bt).sum())
    res[key] = dict(direct_end=d_end, direct_smooth=d_sm, direct_smooth_bias=d_smb, B_P0=B, A1_implied=B / d_smb, mean_lam2=float((lam1 ** 2).mean()))
    print('%s N=%d direct P=%d: end %.5f log-Cesaro %.5f +bias(nu=2) %.5f | B(1) (p<=%d exact) = %.5f -> L(1,Ad) implied = %.5f | mean |lam_p|^2 = %.4f (%.0fs)'
          % (key, N, P, d_end, d_sm, d_smb, P0, B, B / d_smb, (lam1 ** 2).mean(), time.time() - t0))
    sys.stdout.flush()
json.dump(res, open(OUT + 'r6c_g2.json', 'w'), indent=1)
