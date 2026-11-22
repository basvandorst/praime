import math, json
import numpy as np
from sieve import base_primes, primes_in
from sympy import primerange
L0, L1, SEG = 10 ** 9, 2 * 10 ** 9, 10 ** 8
bp = base_primes(int(math.isqrt(L1)) + 2)
tw = []; last = None
for s in range(L0, L1, SEG):
    P = primes_in(s, s + SEG + 2, bp)          # overlap 2 so a twin straddling the edge is seen
    P = P[P < s + SEG + 2]
    t = P[:-1][np.diff(P) == 2]
    t = t[t < s + SEG]
    tw.append(t)
tw = np.concatenate(tw); d = np.diff(tw)
print('twins (lower member in window)', len(tw), 'gaps', len(d), 'all multiples of 6:', bool(np.all(d % 6 == 0)))
v, c = np.unique(d, return_counts=True); N = c.sum(); p = c / N
H = float(-(p * np.log2(p)).sum()); Hmm = H + (len(c) - 1) / (2 * N * math.log(2))
print('H plug-in %.5f  MM %.5f  distinct d %d' % (H, Hmm, len(c)))
# K_tw own: local weights at p>=5 over residues d mod p
K = 0.0
ps = np.array(list(primerange(5, 3 * 10 ** 6)), float)
f = {0: 1 - 2 / ps, 1: 1 - 3 / ps, 2: 1 - 4 / ps}
mean = (f[0] + 2 * f[1] + (ps - 3) * f[2]) / ps
print('mean check vs p(p-2)/(p-1)^2 * (1-2/p)? :', np.allclose(mean, (1 - 2 / ps) ** 2 / (1 - 2 / ps) * 1, rtol=1e-12) if False else (mean[:3], ((ps - 2) ** 2 / ps ** 2)[:3]))
terms = sum(n * (f[i] / mean) * np.log(f[i] / mean) for i, n in ((0, 1), (1, 2), (2, ps - 3))) / ps
K = float(terms.sum()); print('K_tw %.6f nats (last term*p^2 %.3f -> tail ~ %.1e)' % (K, terms[-1] * ps[-1] ** 2, terms[-1] * ps[-1] ** 2 / (ps[-1] * math.log(ps[-1]))))
allp = np.array(list(primerange(3, 10 ** 7)), float); C2 = math.exp(np.log1p(-1 / (allp - 1) ** 2).sum())
for lam in (21.106022426259358,):
    rho = 2 * C2 / lam ** 2
    Ha = (math.log(math.e / (6 * rho)) - K) / math.log(2)
    dd = np.arange(6, 30001, 6, dtype=np.int64); S = np.ones(len(dd))
    for q in primerange(5, 30001):
        r = dd % q; nu = np.where(r == 0, 2, np.where((r == 2) | (r == q - 2), 3, 4))
        S *= (1 - nu / q) / ((q - 2) ** 2 / q ** 2 * q / (q - 2) * (1 - 2 / q)) if False else (1 - nu / q)
    Pd = S * np.exp(-rho * dd); Pd /= Pd.sum(); Hd = float(-(Pd * np.log2(Pd)).sum())
    print('lam %.4f C2 %.8f asym %.5f disc %.5f | data-asym %+.4f data-disc %+.4f' % (lam, C2, Ha, Hd, Hmm - Ha, Hmm - Hd))
# lambda weighted by twin density over the window
x = np.linspace(L0, L1, 10001); wt = 1 / np.log(x) ** 2
print('twin-weighted mean ln x %.4f' % (np.sum(np.log(x) * wt) / wt.sum()))
