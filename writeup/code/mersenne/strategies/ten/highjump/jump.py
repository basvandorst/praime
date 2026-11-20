"""Record gaps ('high jump'): the largest gap below x for the primes (OEIS A002386/A005250 b-files) and for fake primes
up to 10^9 (Cramer: odd numbers with chance 2 / ln n; wheels: numbers without prime factors in W with chance
prod p/(p-1) / ln n; true primes up to max(W, 100) copied).  Predictions: Cramer's model's largest gap below x is about
ln x * ln pi(x) (the largest of pi(x) exponential gaps of mean ln x); Wolf's formula (x / pi(x)) (2 ln pi(x) - ln x + c0),
c0 = ln C2 = -0.4154 (C2 twin constant 0.6601...); Cramer's conjecture ln^2 x.
Prediction written before the run: the primes' largest gaps lie below the Cramer-model fakes and close to Wolf's formula;
fakes with a wheel come closer to the primes.   usage: jump.py out.txt"""
import sys, math
import numpy as np
out = open(sys.argv[1], 'w')
def pr(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')
P = [int(l.split()[1]) for l in open('oeis/b002386.txt') if l.strip() and not l.startswith('#')]
G = [int(l.split()[1]) for l in open('oeis/b005250.txt') if l.strip() and not l.startswith('#')]
from sympy import primepi
C2 = 0.6601618158468696
def wolf(x, pix): return x / pix * (2 * math.log(pix) - math.log(x) + math.log(C2))
pr('records from the OEIS b-files: %d, last p = %d, gap %d; the largest g / ln^2 p: %.4f at p = %d' % (len(P), P[-1], G[-1], max(g / math.log(p) ** 2 for p, g in zip(P, G)), max(zip([g / math.log(p) ** 2 for p, g in zip(P, G)], P))[1]))
xs = [10 ** e for e in range(3, 10)]
pr('     x   primes: largest gap   ln^2 x   ln x ln pi(x)   Wolf')
gm = []
for x in xs:
    g = max(gg for p, gg in zip(P, G) if p + gg <= x); pix = int(primepi(x)); gm.append(g)
    pr('%8.0e %20d %8.1f %15.1f %7.1f' % (x, g, math.log(x) ** 2, math.log(x) * math.log(pix), wolf(x, pix)))
sv_small = np.ones(1001, bool); sv_small[:2] = False
for i in range(2, 32):
    if sv_small[i]: sv_small[i * i::i] = False
small = np.nonzero(sv_small)[0]
rng = np.random.default_rng(23)
def run(W, Xmax=10 ** 9, seg=10 ** 7):
    W = np.array(W); c = float(np.prod(W / (W - 1.0))); cut = max(int(W.max()), 100)
    best = {x: 0 for x in xs}; last = None; gmax = 0
    first = small[small <= cut]
    last = int(first[-1])
    for lo in range(cut + 1, Xmax + 1, seg):
        hi = min(Xmax + 1, lo + seg); n = np.arange(lo, hi)
        ok = np.ones(hi - lo, bool)
        for p in W: ok &= (n % p != 0)
        F = n[ok & (rng.random(hi - lo) < np.minimum(1, c / np.log(n)))]
        if len(F) == 0: continue
        g = np.diff(np.concatenate([[last], F])); ends = F
        # running max over ends, record by x
        cm = np.maximum.accumulate(g); cm = np.maximum(cm, gmax)
        for x in xs:
            if lo <= x < hi:
                k = np.searchsorted(ends, x, 'right') - 1
                best[x] = int(cm[k]) if k >= 0 else gmax
        gmax = int(cm[-1]); last = int(F[-1])
    return [best[x] for x in xs]
for name, W in (('Cramer (odd)', [2]), ('wheel 2..13', [2, 3, 5, 7, 11, 13]), ('wheel primes<=100', list(small[small <= 100]))):
    R = np.array([run(W) for _ in range(6)])
    pr('%-18s largest gap below x = 10^3..10^9, mean of 6 runs: %s' % (name, ' '.join('%.0f' % v for v in R.mean(0))))
    pr('%-18s  spread:                                           %s' % ('', ' '.join('%.0f' % v for v in R.std(0))))
pr('%-18s %s' % ('primes', ' '.join(str(v) for v in gm)))
