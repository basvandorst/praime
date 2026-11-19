"""Width of the Goldbach comet with and without residue-class offsets.  Width: sd over even n in [X/2, X] of r(n)/E(n)
after division by a least-squares quadratic (E = Hardy-Littlewood for the primes, the exact model expectation for
fakes; fakes as in subpoisson3.py).  A random set fills the classes mod small q unevenly, which shifts r(n) for whole
classes of n; the primes fill them evenly.  Offsets removed: for each odd prime q <= Qc, subtract from the residual the
mean over each class n mod q (two passes).  Prediction (from an independent review, written before this run): after
the removal the sieved fakes come much closer to the primes, odd fakes do not change.   usage: subpoisson4.py X runs out.txt"""
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
ns = np.arange(X // 2 + (X // 2) % 2, X + 1, 2)
def I(n):
    v = np.linspace(np.log(2), np.log(n / 2), 4001); t = np.exp(v); f = t / (np.log(t) * np.log(n - t)); return 2 * float(np.sum((f[1:] + f[:-1]) / 2) * (v[1] - v[0]))
C2 = 0.6601618158468696; S = np.ones(X + 1)
for p in P[1:]: S[p::p] *= (p - 1) / (p - 2)
grid = np.linspace(ns[0], ns[-1], 200); Ig = np.array([I(n) for n in grid]); In = np.interp(ns, grid, Ig)
def widths(z):
    z = z / np.polyval(np.polyfit(ns, z, 2), ns) - 1; res = [float(z.std())]
    for Qc in (13, 97):
        y = z.copy()
        for _ in range(2):
            for q in P[(P > 2) & (P <= Qc)]:
                cls = ns % q; mean = np.bincount(cls, y, q) / np.maximum(1, np.bincount(cls, None, q)); y -= mean[cls]
        res.append(float(y.std()))
    return res
r = np.rint(conv(sv.astype(float))); wp = widths(r[ns] / (2 * C2 * S[ns] * In))
pr('X = %d, all %d even n in [X/2, X]; width raw / without class offsets mod odd q <= 13 / <= 97' % (X, len(ns)))
pr('primes                         : %.5f  %.5f  %.5f' % tuple(wp))
rng = np.random.default_rng(44); m = np.arange(X + 1, dtype=float)
for B in (2, 13, 31, 100):
    W = P[P <= B]; ok = np.ones(X + 1, bool); cut = max(B, 100); ok[:cut + 1] = False
    for p in W: ok[::p] = False
    c = float(np.prod(W / (W - 1.0))); pb = np.zeros(X + 1); pb[ok] = np.minimum(1, c / np.log(m[ok])); pb[P[P <= cut]] = 1
    E = conv(pb)[ns]; ws = np.array([widths(np.rint(conv((rng.random(X + 1) < pb).astype(float)))[ns] / E) for _ in range(RUNS)])
    mu, sd = ws.mean(0), ws.std(0)
    pr('fakes avoiding primes up to %3d: %.5f  %.5f  %.5f   (spread %.5f %.5f %.5f, %d runs); primes / fakes %.3f %.3f %.3f' % (B, *mu, *sd, RUNS, *(np.array(wp) / mu)))
