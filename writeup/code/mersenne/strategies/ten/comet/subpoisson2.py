"""Do Goldbach counts fluctuate less than chance?  r(n) = ordered pairs of primes (p, q), p + q = n, even n in [X/2, X].
Reference spread from independent random models: for fakes with chances p_m (true primes up to 100, then numbers
avoiding the multiples of W with chance prod p/(p-1) / ln m), r = 2 sum_{m<n/2} X_m X_{n-m} (+ middle term), so
Var r = 4 sum_{m<n/2} q_m (1 - q_m) with q_m = p_m p_{n-m}, i.e. Var = 2 (conv(p,p) - conv(p^2,p^2)) away from the middle.
Statistic: sd over n of r/E after removing a quadratic drift, divided by the mean of sqrt(Var)/E (the model's own spread).
Prediction written before the run: fakes give 1 (a check of the method); the primes, compared with the spread of the
best model (W up to 13), give clearly less than 1.   usage: subpoisson2.py X runs out.txt"""
import sys
import numpy as np
X, RUNS = int(sys.argv[1]), int(sys.argv[2]); out = open(sys.argv[3], 'w')
def pr(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')
sv = np.ones(X + 1, bool); sv[:2] = False
for i in range(2, int(X ** .5) + 1):
    if sv[i]: sv[i * i::i] = False
P = np.nonzero(sv)[0]; Lf = 1 << (2 * X + 2).bit_length()
def conv(a): F = np.fft.rfft(a, Lf); return np.fft.irfft(F * F, Lf)[:X + 1]
rng = np.random.default_rng(21)
ns = np.sort(rng.choice(np.arange(X // 2 + (X // 2) % 2, X + 1, 2), min(4000, X // 4), replace=False))
def detr(v): return v / np.polyval(np.polyfit(ns, v, 2), ns)
C2 = 0.6601618158468696; S = np.ones(X + 1)
for p in P[1:]: S[p::p] *= (p - 1) / (p - 2)
def I(n):
    t = np.linspace(2, n - 2, 4001); f = 1 / (np.log(t) * np.log(n - t)); return float(np.sum((f[1:] + f[:-1]) / 2) * (t[1] - t[0]))
r = np.rint(conv(sv.astype(float))).astype(np.int64)
HL = np.array([2 * C2 * S[n] * I(n) for n in ns]); qp = detr(r[ns] / HL)
pr('X = %d, %d even n in [X/2, X]; sd of the primes r/HL (drift removed): %.5f; Poisson-for-pairs spread mean sqrt(2/HL): %.5f (ratio %.3f)' % (X, len(ns), qp.std(), np.mean(np.sqrt(2 / HL)), qp.std() / np.mean(np.sqrt(2 / HL))))
m = np.arange(X + 1, dtype=float)
for W in ((2,), (2, 3), (2, 3, 5), (2, 3, 5, 7), (2, 3, 5, 7, 11, 13)):
    ok = np.ones(X + 1, bool); ok[:101] = False
    for p in W: ok[::p] = False
    c = 1.0
    for p in W: c *= p / (p - 1)
    pb = np.zeros(X + 1); pb[ok] = np.minimum(1, c / np.log(m[ok])); pb[P[P <= 100]] = 1
    E = conv(pb)[ns]; V = 2 * (conv(pb) - conv(pb * pb))[ns]; model_sd = np.mean(np.sqrt(V) / E)
    rat = []
    for s in range(RUNS):
        fk = (rng.random(X + 1) < pb).astype(float); rf = np.rint(conv(fk))[ns]
        rat.append(detr(rf / E).std() / model_sd)
    pr('model avoiding %-20s: its own spread %.5f; fakes measured / own spread %.3f (%d runs, sd %.3f); primes / this spread %.3f' % (W, model_sd, np.mean(rat), RUNS, np.std(rat), qp.std() / model_sd))
