import sys, time
import numpy as np
from math import log, log2, sqrt, exp, pi
sys.path.insert(0, "/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins")
from explib import sieve, segments, zeros, primes_upto, spectrum, GAMMA

t0 = time.time()
X = 10**9
min_sqrt = (10**9, 0)
recip = 0.0
logprod = 0.0
prev_last = None
checks = {}
for lo, P in segments(X, 10**8):
    idx = np.flatnonzero(P[:min(10**8, X - lo)]) + lo
    pf = idx.astype(np.float64)
    recip += (1 / pf).sum(); logprod += np.log1p(-1 / pf).sum()
    cum = np.cumsum(P[:min(10**8, X - lo)])
    xs = np.arange(max(lo, 100), lo + len(cum) - 70000, 1000)
    w = np.sqrt(xs).astype(np.int64)
    c = cum[xs + w - lo] - cum[xs - lo]
    i = c.argmin()
    if c[i] < min_sqrt[0]:
        min_sqrt = (int(c[i]), int(xs[i]))
    for d in (10**8, 3 * 10**8, 10**9):
        if lo < d <= lo + len(cum):
            checks[d] = (recip - log(log(d)), logprod + log(log(d)) + GAMMA)
    print(f"  segment {lo:.0e}  [{time.time()-t0:.0f} s]", flush=True)
print(f"E-C19 primes in [x, x + sqrt x] to 1e9 (sampled every 1000): minimum count {min_sqrt[0]} at x = {min_sqrt[1]} (Legendre-type; the conjecture says never 0)")
print("E-C20 Mertens' constants from the primes to 1e9: sum 1/p - ln ln x -> M = 0.2614972, and ln prod(1-1/p) + ln ln x + gamma -> 0:")
for d, (m, e) in checks.items():
    print(f"  x={d:.0e}: M estimate {m:.7f} (error {m - 0.2614972128:+.1e} ~ 1/(x^(1/2) ...)),  Mertens product check {e:+.2e}")

print("E-F8 bits per twin gap against scale (twin gaps in [x, 2x] from the 1e8 sieve):")
isp = sieve(10**8 + 2)
pr = np.flatnonzero(isp[:10**8 + 1]); tw = pr[isp[pr + 2]]
def H(v):
    _, c = np.unique(v, return_counts=True); p = c / c.sum(); return float(-(p * np.log2(p)).sum())
for x in (10**5, 10**6, 10**7, 5 * 10**7):
    a, b = np.searchsorted(tw, [x, 2 * x])
    g = np.diff(tw[a:b]) // 6
    print(f"  x={x:9d}: {len(g)} gaps, H = {H(np.minimum(g, 600)):.3f} bits, log2 of the mean gap {log2(g.mean()*6):.3f}; H(gap | p mod 5) = {H((tw[a:b-1] % 5) * 1000 + np.minimum(g, 600)) - H(tw[a:b-1] % 5):.3f}")

print("E-D13 primes of the form n^2 + 1 for n <= 3e7 against Hardy-Littlewood's (C/2) int dn/ln(n^2+1), C = 1.3728134628:")
N = 3 * 10**7
val = np.ones(N + 1, dtype=bool)
val[0] = False
for q in primes_upto(N):
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
