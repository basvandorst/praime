import sys, time
import numpy as np
import mpmath as mp
from math import log, sqrt, pi, exp
sys.path.insert(0, "/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins")
from explib import zeros, primes_upto, sieve, GAMMA, C2, unfold

t0 = time.time()
gam = zeros()
mp.mp.dps = 20
print("E-X8b Weil's explicit formula, both sides for Gaussian test functions of width s = 3, 8, 20:")
print("  sum_rho h(gamma) = h(i/2) + h(-i/2) + (1/2pi) int h(r) [Re psi(1/4 + ir/2) - ln pi] dr - 2 sum_n Lambda(n) n^(-1/2) g(ln n),  g(u) = (1/2pi) int h(r) e^(-iur) dr")
ps = primes_upto(10**6)
for s_ in (3.0, 8.0, 20.0):
    h = lambda t, s_=s_: mp.exp(-t * t / (2 * s_ * s_))
    g = lambda u, s_=s_: s_ / mp.sqrt(2 * mp.pi) * mp.exp(-s_ * s_ * u * u / 2)
    lhs = 2 * float(sum(mp.exp(-mp.mpf(gm) ** 2 / (2 * s_ * s_)) for gm in gam[gam < 12 * s_]))
    pole = 2 * float(mp.exp(mp.mpf(1) / (8 * s_ * s_)))
    prime_sum = 0.0
    for p in ps:
        q = p
        while q <= 10**6:
            val = float(g(log(q)))
            if val < 1e-30: break
            prime_sum += 2 * log(p) / sqrt(q) * val
            q *= p
    arch = float(mp.quad(lambda t: h(t) * (mp.re(mp.digamma(mp.mpf('0.25') + 1j * t / 2)) - mp.log(mp.pi)), [-15 * s_, -2 * s_, 0, 2 * s_, 15 * s_]) / (2 * mp.pi))
    rhs = pole + arch - prime_sum
    print(f"  s={s_:4.1f}: zeros {lhs:12.6f} | poles {pole:.6f} + archimedean {arch:+.6f} - primes {prime_sum:.6f} = {rhs:12.6f} | difference {lhs - rhs:+.2e}  [{time.time()-t0:.0f} s]", flush=True)

print("E-C18 Lagarias' criterion (RH <=> sigma(n) <= H_n + exp(H_n) ln H_n for all n) on the colossally abundant numbers up to ln N = 2e5:")
events = []
for q in primes_upto(200000).astype(np.float64):
    a_ = 1
    while True:
        eps = log((q ** (a_ + 1) - 1) / (q ** a_ - 1)) / log(q) - 1
        if eps < 1e-7 and a_ > 1: break
        events.append((eps, q)); a_ += 1
        if q ** (a_ + 1) > 1e300: break
events.sort(key=lambda e: -e[0])
lnN = lnsig = 0.0; expo = {}; worst = (-1e9, 0)
for eps, q in events:
    a_ = expo.get(q, 0) + 1; expo[q] = a_
    lnN += log(q); lnsig += log(q ** (a_ + 1) - 1) - (log(q ** a_ - 1) if a_ > 1 else log(q - 1))
    if lnN < 50: continue
    Hn = lnN + GAMMA
    rhs_log = Hn + log(log(Hn))
    margin = (rhs_log - lnsig) * sqrt(lnN) * log(lnN)
    if margin > worst[0] or worst[1] == 0:
        pass
    worst = (margin, lnN) if (worst[1] == 0 or margin < worst[0]) else worst
print(f"  ln(e^H_n ln H_n) - ln sigma(N) > 0 on every CA number; the smallest value of (ln rhs - ln sigma) sqrt(ln N) ln ln N is {worst[0]:.3f} at ln N = {worst[1]:.0f} (Robin's analogue gave 1.48)")

print("E-D12 twin and Sophie Germain count fluctuations in x-space against the zeros (Lambda weights, 1e5 zeros, x in [1e4, 1e6]): E(x) = sum Lambda(n)Lambda(n+h) - S(h) x against -2 S(h) Re sum x^rho/rho")
X = 2 * 10**6 + 2
isp = sieve(X)
pr = np.flatnonzero(isp[:10**6 + 1]); pr = pr[pr > 2]
rho = 0.5 + 1j * gam
xs = np.geomspace(10**4, 10**6, 150)
Zs = np.array([-2 * (2 * (x ** rho / rho).real.sum()) for x in xs])
for name, mask, S in (("twins", isp[pr + 2], 2 * C2), ("Sophie Germain", isp[2 * pr + 1], 2 * C2)):
    v = pr[mask].astype(np.float64)
    w = np.log(v) * (np.log(v + 2) if name == "twins" else np.log(2 * v + 1))
    cum = np.cumsum(w)
    E = np.array([cum[np.searchsorted(v, x, side='right') - 1] for x in xs]) - S * xs
    print(f"  {name:15s}: corr(E, -2 S sum x^rho/rho) = {np.corrcoef(E, S * Zs)[0,1]:+.3f}, regression slope {np.polyfit(S * Zs, E, 1)[0]:.2f} (circle method 1; shot noise sd ~ ln^2 x sqrt(count) = {log(10**6)**2 * sqrt(8169):.0f} against the zeros term ~ {abs(S*Zs[-1]):.0f} at 1e6)")

print("E-B12 Montgomery's F(alpha) from zeros 20000..100000, w(u) = 4/(4 + u^2):")
g = gam[20000:]
T = g[-1]; L = log(T / (2 * pi))
pairs_d = []
for k in range(1, 60):
    pairs_d.append(g[k:] - g[:-k])
dd = np.concatenate(pairs_d)
wgt = 4 / (4 + (dd * L / (2 * pi) * 2 * pi / L) ** 2)
norm = (len(g) * L / (2 * pi))
for alpha in (0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5):
    F = (1 + 2 * (np.cos(alpha * L * dd) * 4 / (4 + dd ** 2)).sum() / len(g)) / 1.0
    print(f"  alpha={alpha:4.2f}: F = {F:.3f}   (Montgomery: {min(alpha, 1.0):.2f} + T^(-2 alpha) ln T ~ {min(alpha,1.0):.2f})")

print("E-X7b Keating-Snaith moments at two heights (|zeta(1/2+it)|^2 and ^4 over 500 units, step 0.5):")
for T0 in (10**4, 5 * 10**4):
    tt = np.arange(T0, T0 + 500, 0.5)
    az = np.array([abs(complex(mp.zeta(mp.mpc(0.5, t)))) for t in tt])
    print(f"  T={T0}: <|zeta|^2> = {(az**2).mean():.3f} (ln(T/2pi) + 2gamma - 1 = {log(T0/(2*pi)) + 2*GAMMA - 1:.3f});  <|zeta|^4> = {(az**4).mean():.1f} ((1/2pi^2) ln^4 T = {log(T0)**4/(2*pi**2):.1f})  [{time.time()-t0:.0f} s]", flush=True)
