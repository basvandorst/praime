"""Width of the Goldbach comet: sd over even n in [X/2, X] of r(n)/E(n) (quadratic drift removed), for the primes
(E = Hardy-Littlewood) and for fakes (E = the exact expectation of the model), the fakes avoiding the multiples of all
primes up to a bound B (true primes up to max(B,100) copied).  The ensemble variance of a model includes the
fluctuation of the total number of fakes, which is common to all n and does not show in the spread over n; so the fair
reference is the measured spread of the fakes, not the formula (checked in subpoisson2.py: 0.8 to 0.94 of it).
Prediction written before the run: the ratio primes / fakes rises towards 1 as B grows.  usage: subpoisson3.py X runs out.txt"""
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
rng = np.random.default_rng(33)
ns = np.sort(rng.choice(np.arange(X // 2 + (X // 2) % 2, X + 1, 2), min(4000, X // 4), replace=False))
def detr(v): return v / np.polyval(np.polyfit(ns, v, 2), ns)
C2 = 0.6601618158468696; S = np.ones(X + 1)
for p in P[1:]: S[p::p] *= (p - 1) / (p - 2)
def I(n):
    t = np.linspace(2, n - 2, 4001); f = 1 / (np.log(t) * np.log(n - t)); return float(np.sum((f[1:] + f[:-1]) / 2) * (t[1] - t[0]))
r = np.rint(conv(sv.astype(float))).astype(np.int64)
HL = np.array([2 * C2 * S[n] * I(n) for n in ns]); sp = detr(r[ns] / HL).std()
pr('X = %d, %d even n in [X/2, X]; primes: sd of r/HL %.5f' % (X, len(ns), sp))
m = np.arange(X + 1, dtype=float)
for B in (2, 13, 31, 100):
    W = P[P <= B]; ok = np.ones(X + 1, bool); cut = max(B, 100); ok[:cut + 1] = False
    for p in W: ok[::p] = False
    c = float(np.prod(W / (W - 1.0)))
    pb = np.zeros(X + 1); pb[ok] = np.minimum(1, c / np.log(m[ok])); pb[P[P <= cut]] = 1
    E = conv(pb)[ns]; sds = []
    for s in range(RUNS):
        fk = (rng.random(X + 1) < pb).astype(float); sds.append(detr(np.rint(conv(fk))[ns] / E).std())
    sds = np.array(sds)
    pr('fakes avoiding the primes up to %3d (chance factor %.2f, largest chance %.2f): sd %.5f (runs %d, spread %.5f); primes / fakes %.3f' % (B, c, pb[cut + 1:].max(), sds.mean(), RUNS, sds.std(), sp / sds.mean()))
