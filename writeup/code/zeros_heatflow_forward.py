import sys, time
import numpy as np
from math import log, sqrt, pi
from scipy.special import expi
sys.path.insert(0, "/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins")
from explib import zeros, unfold, sieve

t0 = time.time()
gam = zeros(4000)
print("E-K1 forward de Bruijn-Newman flow dx_j/dt = sum_k 2/(x_j - x_k) on zeros 500..3500 (edges 500 each side frozen), unfolded spacing variance vs t (GUE 0.18, lattice 0):")
x = gam.copy()
core = slice(500, 3500)
def var_core(x):
    s = unfold(x[core]); return s.var() / s.mean() ** 2
print(f"  t=0.00: var {var_core(x):.4f}")
t = 0.0
for target in (0.02, 0.05, 0.1, 0.2):
    while t < target:
        diff = x[:, None] - x[None, :]
        np.fill_diagonal(diff, np.inf)
        v = (2 / diff).sum(axis=1)
        v[:500] = 0; v[3500:] = 0
        gaps = np.diff(x)
        step = min(target - t, 0.02 * (gaps[500:3500] ** 2).min() / 4)
        x += step * v
        t += step
    print(f"  t={t:.2f}: var {var_core(x):.4f}  (mean spacing of the core {np.diff(x[core]).mean():.4f} vs {np.diff(gam[core]).mean():.4f} at t = 0)  [{time.time()-t0:.0f} s]", flush=True)

print("E-C4b Schoenfeld's pi-li inequality restricted to its range x >= 2657, up to 1e7:")
X = 10**7
isp = sieve(X)
pic = np.cumsum(isp)
xs = np.arange(2657, X + 1, dtype=np.float64)
li = expi(np.log(xs)) - expi(log(2))
r = np.abs(pic[2657:] - li) / (np.sqrt(xs) * np.log(xs) / (8 * pi))
i = np.argmax(r)
print(f"  max ratio {r[i]:.4f} at x = {int(xs[i])}; ratio at 1e5, 1e6, 1e7: {r[10**5-2657]:.4f}, {r[10**6-2657]:.4f}, {r[10**7-2657]:.4f}")

print("E-F5b digit-sum parity: is the 0.495 for primes in [1e7, 1e8) a mod-3 artefact? Compare all n coprime to 6 in the same range:")
rng = np.random.default_rng(5)
n = rng.integers(10**7, 10**8, size=600000)
n = n[(n % 2 == 1) & (n % 3 != 0)]
s = np.array([sum(map(int, str(v))) for v in n])
print(f"  n coprime to 6: P(even digit sum) = {(s % 2 == 0).mean():.4f} +- {sqrt(0.25/len(n)):.4f} (primes gave 0.4950 +- 0.0007)")
m = n % 5 != 0
print(f"  n coprime to 30: {(s[m] % 2 == 0).mean():.4f} +- {sqrt(0.25/m.sum()):.4f}")
print(f"[{time.time()-t0:.0f} s]")
