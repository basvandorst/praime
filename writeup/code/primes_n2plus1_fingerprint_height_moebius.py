import sys, time
import numpy as np
from math import log, log2, sqrt, exp, pi
sys.path.insert(0, "/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins")
from explib import sieve, segments, zeros, primes_upto, spectrum, GAMMA

t0 = time.time()
print("E-D13 primes of the form n^2 + 1 for n <= 3e7 against Hardy-Littlewood's (C/2) int dn/ln(n^2+1), C = 1.3728134628:")
N = 3 * 10**7
val = np.ones(N + 1, dtype=bool)
val[0] = False
for q in primes_upto(N).tolist():
    if q == 2:
        val[1::2] = val[1::2] & False if False else val[1::2]
        continue
    if q % 4 != 1:
        continue
    r = pow(int(q) - 1, (int(q) + 1) // 4, int(q)) if q % 8 == 5 else None
    rr = None
    for cand in range(2, 100):
        if pow(cand, (q - 1) // 2, q) == q - 1:
            rr = pow(cand, (q - 1) // 4, q); break
    if rr is None:
        continue
    for root in (rr, q - rr):
        start = root
        if start * start + 1 == q:
            start += q
        val[start::q] = False
val[1::2] = False
val[1] = True
cnt = int(val[1:N + 1].sum())
n = np.arange(2, N + 1, dtype=np.float64)
hl = 0.5 * 1.3728134628 * (1 / np.log(n * n + 1)).sum()
print(f"  count {cnt} (conjecture 58 page: 456362 at 1e7), HL {hl:.0f}, ratio {cnt/hl:.4f}  [{time.time()-t0:.0f} s]")

print("E-B13 the fingerprint law against height: mean cos(gamma ln 2) in five bins of 20000 zeros against -Lambda(2)/(sqrt 2 ln(T/2pi e)) with T the bin's mean height:")
gam = zeros()
for k in range(5):
    g = gam[k * 20000:(k + 1) * 20000]
    pred = -log(2) / sqrt(2) / log(g.mean() / (2 * pi * np.e))
    print(f"  zeros {k*20000}..{(k+1)*20000} (T ~ {g.mean():.0f}): measured {np.cos(g * log(2)).mean():+.4f}, predicted {pred:+.4f}")

print("E-A14 blind zeros from the Moebius function alone: peaks of |sum mu(n) n^(-1/2-it)| (n <= 1e6) on [10, 60] step 0.02:")
mu = np.ones(10**6 + 1, dtype=np.int8); mu[0] = 0
for q in primes_upto(10**6):
    mu[q::q] *= -1
    if q <= 1000: mu[q*q::q*q] = 0
nn = np.arange(1, 10**6 + 1, dtype=np.float64)
tg = np.arange(10.0, 60.0, 0.02)
S = np.abs(spectrum(nn, mu[1:].astype(np.float64), tg))
found = []
for i in range(1, len(tg) - 1):
    if S[i] > S[i - 1] and S[i] > S[i + 1] and S[i] > 8:
        a, b, c = S[i - 1], S[i], S[i + 1]
        found.append(tg[i] + 0.5 * (a - c) / (a - 2 * b + c) * 0.02)
found = np.array(found)
true = gam[(gam > 10) & (gam < 60)]
d = np.array([np.abs(found - z).min() for z in true])
print(f"  {len(found)} peaks above 8, zeros within 0.1: {(d < 0.1).sum()}/{len(true)}, rms {sqrt((d[d<0.5]**2).mean()):.4f}, spurious {sum(np.abs(true - f).min() > 0.1 for f in found)}  [{time.time()-t0:.0f} s]")
