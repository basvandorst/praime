import sys, time
import numpy as np
import mpmath as mp
from math import log, sqrt, pi, exp
sys.path.insert(0, "/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins")
from explib import zeros, unfold, primes_upto, GAMMA

t0 = time.time()
gam = zeros()
print("E-B10 prime peaks in the form factor of the zeros: K(tau) = |sum_j exp(2 pi i tau x_j)|^2 / N with unfolded x_j (zeros 20000..100000); peaks expected at tau_p = ln p / (2 pi d), d the mean density ln(T/2pi)/(2pi):")
g = gam[20000:]
x = np.cumsum(np.concatenate([[0.0], unfold(g)]))
d = np.log(g.mean() / (2 * pi)) / (2 * pi)
taus = np.arange(0.02, 0.45, 0.0005)
K = np.array([abs(np.exp(2j * pi * t * x).sum()) ** 2 / len(x) for t in taus])
print(f"  background: median K over the grid {np.median(K):.3f}, median tau {np.median(taus):.3f} (GUE K = tau)")
for p_ in (2, 3, 4, 5, 7, 8, 9, 6, 10):
    tp = log(p_) / (2 * pi * d)
    kp = abs(np.exp(2j * pi * tp * x).sum()) ** 2 / len(x)
    kn = np.mean([abs(np.exp(2j * pi * (tp + dd) * x).sum()) ** 2 / len(x) for dd in (-0.003, -0.002, 0.002, 0.003)])
    print(f"  m={p_:2d}: tau_m = {tp:.4f}, K(tau_m) = {kp:8.2f}, K nearby = {kn:6.2f}, ratio {kp/kn:6.1f}   {'prime power' if p_ in (2,3,4,5,7,8,9) else 'composite'}")
print("  predicted tau_p for p = 2, 3, 4, 5, 7, 8, 9:", ", ".join(f"{log(p)/(2*pi*d):.4f}" for p in (2, 3, 4, 5, 7, 8, 9)), f"(d = {d:.3f})")
print("  GUE background K(tau) = tau; the peak heights follow Lambda(p)^2/p * N/(... ) in Berry's semiclassical picture")

print("E-B11 number-variance plateau against height: Sigma^2(L) for L = 5, 10, 20 on zeros 10000..30000 and 80000..100000:")
for lo, hi in ((10000, 30000), (50000, 70000), (80000, 100000)):
    xx = np.cumsum(np.concatenate([[0.0], unfold(gam[lo:hi])]))
    line = f"  zeros {lo}..{hi} (T ~ {gam[(lo+hi)//2]:.0f}):"
    for Lw in (2, 5, 10, 20):
        starts = np.arange(0, xx[-1] - Lw, Lw / 3)
        c = np.searchsorted(xx, starts + Lw) - np.searchsorted(xx, starts)
        line += f"  L={Lw}: {c.var():.3f}"
    print(line)

print("E-X8 Weil's explicit formula as an identity: for h(t) = exp(-t^2/(2 s^2)) with s = 3, sum_gamma h(gamma) (both signs, 1e5 zeros) against h-hat(0) ln(pi) terms... computed as: sum over zeros = h^(i/2)+h^(-i/2) - sum_n Lambda(n)/sqrt(n) (g(ln n) + g(-ln n)) + (1/2pi) int h(t) Re psi(1/4 + it/2) dt - h-hat ln pi, with g the Fourier transform of h")
mp.mp.dps = 20
s_ = 3.0
h = lambda t: mp.exp(-t * t / (2 * s_ * s_))
gfun = lambda u: s_ / mp.sqrt(2 * mp.pi) * mp.exp(-s_ * s_ * u * u / 2)
lhs = float(sum(mp.exp(-mp.mpf(gm) ** 2 / (2 * s_ * s_)) for gm in gam[:20000])) * 2
pole = float((h(mp.mpc(0, 0.5)) + h(mp.mpc(0, -0.5))).real)
ps = primes_upto(200000)
prime_sum = 0.0
for p in ps:
    q = p
    while q <= 200000:
        prime_sum += float(log(p) / sqrt(q) * 2 * gfun(log(q)))
        q *= p
arch = float(mp.quad(lambda t: h(t) * (mp.re(mp.digamma(mp.mpf('0.25') + 1j * t / 2)) - mp.log(mp.pi)), [-40, 40]) / (2 * mp.pi))
rhs = pole - prime_sum + arch
print(f"  sum over zeros = {lhs:.6f};  pole terms {pole:.6f} - prime sum {prime_sum:.6f} + archimedean {arch:.6f} = {rhs:.6f};  difference {lhs - rhs:+.2e}  [{time.time()-t0:.0f} s]")
