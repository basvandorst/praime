"""Trains of consecutive primes: f(n) = number of ways to write n as a sum of one or more consecutive primes.
Moser (1963): (1/x) sum_{n <= x} f(n) -> ln 2.  Measured: the mean, the mean over even and over odd n, the distribution
of f(n) (n in [x/2, x]) against Poisson with the same mean, and the same for Cramer fakes (odd fakes, true primes up to
100).  Prediction written before the run: the mean is near ln 2 from above; f is close to Poisson.   usage: trains.py x out.txt"""
import sys, math
import numpy as np
x = int(sys.argv[1]); out = open(sys.argv[2], 'w')
def pr(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')
sv = np.ones(x + 1, bool); sv[:2] = False
for i in range(2, int(x ** .5) + 1):
    if sv[i]: sv[i * i::i] = False
def trains(P):
    f = np.zeros(x + 1, np.int32); cs = np.concatenate([[0], np.cumsum(P.astype(np.int64))])
    for L in range(1, len(P) + 1):                    # length L: sums cs[i + L] - cs[i]
        s = cs[L:] - cs[:-L]
        s = s[s <= x]
        if len(s) == 0: break
        np.add.at(f, s, 1)
    return f
P = np.nonzero(sv)[0]
f = trains(P)
n = np.arange(x + 1); lo = x // 2
pr('x = %d: mean f over 1..x %.5f (ln 2 = %.5f); over even n %.5f, over odd n %.5f; over odd composite n %.5f' % (x, f[1:].mean(), math.log(2), f[2::2].mean(), f[1::2].mean(), f[1::2][~sv[1::2]].mean()))
seg = f[lo:]; mu = seg.mean(); K = int(seg.max()) + 1; c = np.bincount(seg, minlength=K)
pr('n in [x/2, x]: mean %.4f; counts of f = 0..%d: %s' % (mu, K - 1, ' '.join(str(v) for v in c)))
pr('Poisson with the same mean:          %s' % ' '.join('%.0f' % (len(seg) * math.exp(-mu) * mu ** k / math.factorial(k)) for k in range(K)))
for name, sel in (('even n', seg[(np.arange(lo, x + 1) % 2) == 0]), ('odd n', seg[(np.arange(lo, x + 1) % 2) == 1])):
    m = sel.mean(); cc = np.bincount(sel, minlength=K)
    pr('%-7s mean %.4f var/mean %.4f; counts %s; Poisson %s' % (name, m, sel.var() / m, ' '.join(str(v) for v in cc), ' '.join('%.0f' % (len(sel) * math.exp(-m) * m ** k / math.factorial(k)) for k in range(K))))
for q in (3, 5, 7):
    pr('mean f by n mod %d (n in [x/2, x]): %s' % (q, ' '.join('%.4f' % seg[(np.arange(lo, x + 1) % q) == r].mean() for r in range(q))))
pr('largest f up to x: %d, first at n = %d' % (int(f.max()), int(np.argmax(f))))
rng = np.random.default_rng(18); m2 = np.arange(x + 1, dtype=float); pb = np.zeros(x + 1); pb[101::2] = np.minimum(1, 2 / np.log(m2[101::2]))
for r in range(2):
    Fk = np.concatenate([P[P <= 100], np.nonzero(rng.random(x + 1) < pb)[0]]); ff = trains(Fk); s2 = ff[lo:]
    pr('odd fakes run %d: mean over 1..x %.5f; [x/2, x] mean %.4f var/mean %.4f; even n mean %.4f, odd n mean %.4f; by n mod 3: %s' % (r, ff[1:].mean(), s2.mean(), s2.var() / s2.mean(), s2[(np.arange(lo, x + 1) % 2) == 0].mean(), s2[(np.arange(lo, x + 1) % 2) == 1].mean(), ' '.join('%.4f' % s2[(np.arange(lo, x + 1) % 3) == rr].mean() for rr in range(3))))
