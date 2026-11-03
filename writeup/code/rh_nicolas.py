import sys, time
import numpy as np
from math import log, sqrt, exp

GAMMA = 0.5772156649015329
X = int(float(sys.argv[2])) if len(sys.argv) > 2 else 10**8
gam = np.loadtxt(sys.argv[1], max_rows=100000)
rho = 0.5 + 1j * gam
t0 = time.time()
s = np.ones(X + 1, dtype=bool); s[:2] = False
for i in range(2, int(X**0.5) + 1):
    if s[i]: s[i * i::i] = False
p = np.flatnonzero(s).astype(np.float64)
lnp = np.log(p)
theta = np.cumsum(lnp)
lnNphi = np.cumsum(-np.log1p(-1.0 / p))
margin = lnNphi - GAMMA - np.log(np.log(theta))
print(f"primes to {X:.0e}: {len(p)}; Nicolas margin m_k = ln(N_k/phi(N_k)) - gamma - ln ln N_k  [{time.time()-t0:.0f} s]")
k0 = 2
m = margin[k0:]
print(f"  min over k >= 3: {m.min():.3e} at p_k = {int(p[k0 + m.argmin()])}; negative values: {(m < 0).sum()}  (RH <=> all positive, Nicolas 1983)")
scaled = margin * np.sqrt(theta) * np.log(p)
print("  scaled margin m_k sqrt(ln N_k) ln p_k at decades (explicit-formula leading term 2, oscillation from the zeros):")
print("     p_k         m_k         scaled   zeros prediction 2 + osc")
for d in (3, 4, 5, 6, 7, 8):
    x = 10**d
    if x > X: break
    k = np.searchsorted(p, x) - 1
    xk = p[k]
    osc = (2 * ((xk ** (rho - 0.5)) / (rho * (1 - rho))).real.sum())
    print(f"  {int(xk):10d}  {margin[k]:.4e}  {scaled[k]:8.4f}   {2 + osc:8.4f}")
sub = np.arange(100, len(p), max(1, len(p) // 4000))
xs = p[sub]
osc = np.array([2 * ((x ** (rho - 0.5)) / (rho * (1 - rho))).real.sum() for x in xs])
sc = scaled[sub]
a, b = np.polyfit(osc, sc, 1)
print(f"  over {len(sub)} sample points p_k in [{int(xs[0])}, {int(xs[-1])}]: scaled margin = {b:.3f} + {a:.3f} * osc(zeros)   (prediction 2 + 1.0 * osc); corr = {np.corrcoef(osc, sc)[0,1]:+.3f}")
print(f"  slack: the constant 2 against an oscillation of amplitude at most sum 1/|rho|^2 = {2*(1/np.abs(rho)**2).sum():.4f} (converges to 0.0462): RH-violating zeros would need x^(beta - 1/2) to reach that constant")
