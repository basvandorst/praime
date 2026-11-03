import sys, time
import numpy as np
from math import log, sqrt, pi, exp
sys.path.insert(0, "/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins")
from explib import zeros, unfold, smooth_N, primes_upto

GAMMA = 0.5772156649015329
gam = zeros()
t0 = time.time()
g = gam[20000:]
s = unfold(g)
print(f"zeros 20000..100000: mean spacing {s.mean():.4f}")

print("E-B1 number variance Sigma^2(L) of the unfolded zeros vs GUE (2/pi^2)(ln(2 pi L) + gamma + 1 - pi^2/8):")
x = np.cumsum(np.concatenate([[0.0], s]))
for Lw in (0.5, 1, 2, 5, 10, 20, 50):
    starts = np.arange(0, x[-1] - Lw, Lw / 3)
    counts = np.searchsorted(x, starts + Lw) - np.searchsorted(x, starts)
    gue = 2 / pi**2 * (log(2 * pi * Lw) + GAMMA + 1 - pi**2 / 8)
    print(f"  L={Lw:5.1f}: var {counts.var():.3f}   GUE {gue:.3f}   Poisson {Lw:.1f}")

print("E-B2 pair correlation R2(r) in bins vs 1 - (sin pi r / pi r)^2:")
n = len(x)
R = []
for i in range(1, 40):
    R.append(x[i:] - x[:-i])
d = np.concatenate(R)
edges = np.arange(0, 3.01, 0.25)
h = np.histogram(d, bins=edges)[0] / (n * 0.25)
for a, b, v in zip(edges[:-1], edges[1:], h):
    m = 0.5 * (a + b)
    gue = 1 - (np.sin(pi * m) / (pi * m)) ** 2
    print(f"  [{a:4.2f},{b:4.2f}) {v:.3f}  GUE {gue:.3f}")

print("E-B5 S(t) = N(t) - smooth(t): mean and variance against Selberg's (1/2pi^2) ln ln t:")
Nt = np.arange(1, len(gam) + 1)
S = Nt - 0.5 - smooth_N(gam)
for lo, hi in ((1000, 10000), (10000, 50000), (50000, 100000)):
    seg = S[lo:hi]
    print(f"  zeros {lo}..{hi}: mean {seg.mean():+.3f}, var {seg.var():.4f}, Selberg {1/(2*pi**2)*log(log(gam[(lo+hi)//2])):.4f}, skew {((seg-seg.mean())**3).mean()/seg.std()**3:+.3f}, kurt {((seg-seg.mean())**4).mean()/seg.std()**4:.3f}")

print("E-B7 level repulsion: fraction of spacings below s0 vs GUE (~ (pi^2/9) s0^3 for small s0) and Poisson (s0):")
for s0 in (0.05, 0.1, 0.2, 0.3):
    print(f"  s0={s0}: zeros {(s < s0).mean():.5f}   GUE small-s {pi**2/9*s0**3:.5f}   Poisson {1-exp(-s0):.5f}")

print("E-B8 spacing autocorrelation of the zeros vs a GUE matrix sample (numpy eigenvalues, 3 matrices of size 3000, bulk):")
def acf(v, k):
    v = v - v.mean(); return (v[:-k] * v[k:]).mean() / v.var()
print("  zeros: lag1 %.4f lag2 %.4f lag3 %.4f" % (acf(s, 1), acf(s, 2), acf(s, 3)))
rng = np.random.default_rng(3)
vals = []
for _ in range(3):
    N = 3000
    A = rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))
    H = (A + A.conj().T) / 2
    e = np.linalg.eigvalsh(H)
    e = e[N // 4: 3 * N // 4]
    dens = np.sqrt(4 * N - (e * np.sqrt(2)) ** 2) / (2 * pi)
    sp = np.diff(e) * np.sqrt(2) * dens[:-1] / np.sqrt(2)
    sp = sp / sp.mean()
    vals.append((acf(sp, 1), acf(sp, 2), acf(sp, 3)))
vals = np.array(vals)
print("  GUE:   lag1 %.4f lag2 %.4f lag3 %.4f (sd %.3f)" % (vals[:, 0].mean(), vals[:, 1].mean(), vals[:, 2].mean(), vals[:, 0].std()))

print("E-H4 arithmetic fingerprint: mean of cos(gamma ln m) over the first 100000 zeros (the explicit formula predicts a negative mean of order -(ln p)/sqrt(p) * N/(2 pi ...) only when m is a prime power):")
for m in (2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13):
    c = np.cos(gam * log(m)).mean()
    print(f"  m={m:2d}: {c:+.4f}   {'prime power' if m in (2,3,4,5,7,8,9,11,13) else 'composite'}")
print(f"[{time.time()-t0:.0f} s]")
