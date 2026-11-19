"""Goldbach comet: r(n) = number of ordered pairs of primes (p, q) with p + q = n, for even n <= X (numpy FFT, rounded,
checked by direct count on a sample). Against Hardy-Littlewood: r(n) ~ 2 C2 S(n) I(n), S(n) = prod_{p | n, p > 2} (p-1)/(p-2),
I(n) = integral_2^{n-2} dt / (ln t ln(n-t)).  Predictions written before the run: the classes 3 | n and 3 not | n differ by
a factor 2; 5: 4/3; 7: 6/5; after division by S(n) the bands merge.   usage: comet.py X out.txt"""
import sys, math
import numpy as np
X = int(sys.argv[1]); out = open(sys.argv[2], 'w')
def pr(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')
sv = np.ones(X + 1, bool); sv[:2] = False
for i in range(2, int(X ** .5) + 1):
    if sv[i]: sv[i * i::i] = False
L = 1 << (2 * X + 2).bit_length()
F = np.fft.rfft(sv.astype(float), L); r = np.rint(np.fft.irfft(F * F, L)[:X + 1]).astype(np.int64)
P = np.nonzero(sv)[0]
rng = np.random.default_rng(0)
for n in list(rng.choice(np.arange(4, X + 1, 2), 30)) + [4, 6, X - X % 2]:
    n = int(n); assert r[n] == int(np.sum(sv[n - P[P < n]])), n
pr('X = %d; r(n) checked against a direct count at 33 even n' % X)
E = np.arange(6, X + 1, 2)
pr('every even n from 4 to X has r(n) >= 1:', bool(np.all(r[np.arange(4, X + 1, 2)] >= 1)))
C2 = 1.0
for p in P[1:]:
    if p > 10 ** 6: break
    C2 *= 1 - 1 / (p - 1) ** 2
C2 = 0.6601618158468696                       # twin prime constant (literature); the partial product above:
pr('twin prime constant C2 = 0.66016181584687 (literature); partial product over odd p <= min(X,10^6): %.8f' % (np.prod([1 - 1 / (p - 1) ** 2 for p in P[1:] if p <= 10 ** 6])))
oddp = P[1:]
S = np.ones(X + 1)
for p in oddp:
    S[p::p] *= (p - 1) / (p - 2)
def I(n):
    t = np.linspace(2, n - 2, 4001); f = 1 / (np.log(t) * np.log(n - t)); return float(np.sum((f[1:] + f[:-1]) / 2) * (t[1] - t[0]))
pr('\n1. Bands. Mean of r(n) / (n / ln^2 n) over even n in (X/2, X], by class; ratio of the classes; Hardy-Littlewood ratio')
half = E[E > X // 2]; base = r[half] / (half / np.log(half) ** 2)
for p in (3, 5, 7, 11, 13):
    a = base[half % p == 0].mean(); b = base[half % p != 0].mean()
    pr('  p = %2d: divisible %.4f, not %.4f, ratio %.4f; predicted %.4f' % (p, a, b, a / b, (p - 1) / (p - 2)))
pr('\n2. Ratio r(n) / (2 C2 S(n) I(n)) at even n (trapezoid with 4001 points)')
qs = [n for n in (100, 1000, 5000, 10 ** 4, 10 ** 5, 10 ** 6, 2 * 10 ** 6) if n <= X]
for n in qs:
    pr('  n = %8d: r = %6d, prediction %9.2f, ratio %.4f' % (n, r[n], 2 * C2 * S[n] * I(n), r[n] / (2 * C2 * S[n] * I(n))))
pr('\n3. Spread of the ratio over windows of even n (mean, standard deviation; Poisson-like spread 1/sqrt(r) for comparison)')
for lo, hi in ((1000, 5000), (X // 2, X)):
    ns = np.arange(lo + (lo % 2), hi + 1, 2)
    if len(ns) > 4000: ns = np.sort(rng.choice(ns, 4000, replace=False))
    q = np.array([r[n] / (2 * C2 * S[n] * I(n)) for n in ns]); rr = r[ns]
    pr('  even n in [%d, %d] (%d values): mean %.4f, sd %.4f; mean 1/sqrt(r) %.4f; within the bands (raw r / (n/ln^2 n)) sd %.4f' % (lo, hi, len(ns), q.mean(), q.std(), np.mean(1 / np.sqrt(rr)), np.std(rr / (ns / np.log(ns) ** 2))))
pr('\n4. The lowest points: the even n in (X/2, X] with the smallest r(n), and S(n)')
idx = half[np.argsort(r[half])[:5]]
for n in idx: pr('  n = %d = %s, r = %d, S(n) = %.3f' % (n, 'x'.join(str(int(p)) for p in P[:30] if n % p == 0), r[n], S[n]))
