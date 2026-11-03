import sys, time
import numpy as np
import mpmath as mp
from math import log, sqrt, pi, exp
from scipy.special import expi
sys.path.insert(0, "/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins")
from explib import sieve, zeros, lam_array, C2, GAMMA

t0 = time.time()
LZ = np.load("/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins/lzeros.npy", allow_pickle=True)
g3, g4 = LZ[0], LZ[1]
X = 10**8
isp = sieve(X + 8)
pr = np.flatnonzero(isp[:X + 1])
lnp = np.log(pr[1:].astype(np.float64))

print("E-R1 Chebyshev bias mod 4 from the 170 zeros of L(s, chi4): theta(x;4,3) - theta(x;4,1) = sum_rho x^rho/rho (over chi4 zeros, both signs) + lower order; test at x in [1e4, 1e8]")
p = pr[1:]
d4 = np.cumsum(np.where(p % 4 == 3, lnp, -lnp))
xs = np.geomspace(10**4, X, 400)
idx = np.searchsorted(p, xs, side='right') - 1
real4 = d4[idx] / np.sqrt(xs)
r4 = 0.5 + 1j * g4
pred4 = np.array([2 * ((x ** r4) / r4).real.sum() for x in xs]) / np.sqrt(xs)
print(f"  corr(real, zeros) = {np.corrcoef(real4, pred4)[0,1]:+.3f}, slope {np.polyfit(pred4, real4, 1)[0]:.3f}, mean real {real4.mean():+.3f}, mean pred {pred4.mean():+.3f}  (a residual constant is the prime-square term: theta counts only primes, the squares of primes are all 1 mod 4)")
print("E-R2 the same mod 3 with the 158 zeros of L(s, chi3): theta(x;3,2) - theta(x;3,1)")
d3 = np.cumsum(np.where(p % 3 == 2, lnp, -lnp))
real3 = d3[idx] / np.sqrt(xs)
r3 = 0.5 + 1j * g3
pred3 = np.array([2 * ((x ** r3) / r3).real.sum() for x in xs]) / np.sqrt(xs)
print(f"  corr(real, zeros) = {np.corrcoef(real3, pred3)[0,1]:+.3f}, slope {np.polyfit(pred3, real3, 1)[0]:.3f}, mean real {real3.mean():+.3f}, mean pred {pred3.mean():+.3f}  [{time.time()-t0:.0f} s]")

print("E-R3 twin races: lesser twin members by class mod 5 (classes 1, 2, 4 of p mod 5) and mod 4 (1, 3), to 1e8; share of x where each class leads, and the deviation from equal shares in units of sqrt(count)")
tw = pr[isp[pr + 2]]; tw = tw[tw > 5]
for mod, classes in ((5, (1, 2, 4)), (4, (1, 3))):
    counts = np.vstack([np.cumsum(tw % mod == c) for c in classes]).astype(np.float64)
    mean = counts.sum(axis=0) / len(classes)
    dev = (counts - mean) / np.sqrt(np.maximum(mean, 1))
    lead = np.argmax(counts, axis=0)
    print(f"  mod {mod}: final counts {counts[:, -1].astype(int).tolist()}; leads: " + ", ".join(f"class {c}: {100*(lead==i).mean():.1f} %" for i, c in enumerate(classes)) + "; max |dev| " + ", ".join(f"{np.abs(dev[i]).max():.2f}" for i in range(len(classes))) + "; mean dev " + ", ".join(f"{dev[i, 1000:].mean():+.2f}" for i in range(len(classes))))

print("E-A11 Riemann's 1859 formula pi(x) = R(x) - sum_rho R(x^rho) (R the Riemann function), 1000 zeros, x in [1e3, 1e6]:")
mp.mp.dps = 15
gam = zeros(1000)
rho = 0.5 + 1j * gam
pic = np.cumsum(isp[:10**6 + 1])
xs2 = np.geomspace(1000, 10**6, 200)
def R(z):
    return mp.riemannr(z)
realpi = np.array([pic[int(x)] for x in xs2])
Rx = np.array([float(R(x)) for x in xs2])
osc = np.array([-2 * sum((R(mp.mpf(x) ** mp.mpc(0.5, g))).real for g in gam[:300]) for x in xs2[::4]])
print(f"  mean |pi - R| = {np.abs(realpi - Rx).mean():.2f}; with 300 zeros: mean |pi - R + sum R(x^rho)| = {np.abs(realpi[::4] - Rx[::4] - osc).mean():.2f}; corr(pi - R, zeros term) = {np.corrcoef(realpi[::4] - Rx[::4], osc)[0,1]:+.3f}  [{time.time()-t0:.0f} s]")

print("E-D6 prime quadruplets (p, p+2, p+6, p+8) to 1e8 against Hardy-Littlewood with the exact constant:")
quad = pr[isp[pr + 2] & isp[pr + 6] & isp[pr + 8]]
c4 = 1.0
for q in np.flatnonzero(sieve(10**6)):
    nu = len(set(o % q for o in (0, 2, 6, 8)))
    c4 *= (1 - nu / q) / (1 - 1 / q) ** 4
t = np.linspace(11, X, 200001)
hl = c4 * np.trapezoid(1 / np.log(t) ** 4, t)
print(f"  {len(quad)} quadruplets (first {quad[:4].tolist()}); HL constant {c4:.5f} (literature 4.15118); HL integral {hl:.0f}; ratio {len(quad)/hl:.4f}")

print("E-D7 Brun's constant: sum of 1/p + 1/(p+2) over twins to 1e8, plus the Hardy-Littlewood tail 2C2 * 2 / ln x:")
tw_all = pr[isp[pr + 2]].astype(np.float64)
B = (1 / tw_all + 1 / (tw_all + 2)).sum()
tail = 4 * C2 / log(X)
print(f"  partial sum {B:.6f}, tail estimate {tail:.6f}, B2 ~ {B + tail:.6f} (Nicely: 1.902160583)  [{time.time()-t0:.0f} s]")
