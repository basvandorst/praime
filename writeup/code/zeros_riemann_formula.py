import sys, time
import numpy as np
import mpmath as mp
from math import log
from scipy.special import expi
sys.path.insert(0, "/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins")
from explib import sieve, zeros

t0 = time.time()
mp.mp.dps = 15
gam = zeros(300)
isp = sieve(10**6 + 1)
pic = np.cumsum(isp)
xs = np.geomspace(1000, 10**6, 120)
realpi = np.array([pic[int(x)] for x in xs])
li = expi(np.log(xs))
print("E-A11 Riemann's explicit formula for pi(x): pi(x) = R(x) - sum_rho li(x^rho) + ..., R(x) = sum mu(n)/n li(x^(1/n)) (mpmath riemannr), li(x^rho) = Ei(rho ln x); 300 zeros, x in [1e3, 1e6]")
Rx = np.array([float(mp.riemannr(x)) for x in xs])
zsum = np.array([float(-2 * sum((mp.ei(mp.mpc(0.5, g) * log(x))).real for g in gam)) for x in xs])
pred = Rx + zsum
print(f"  mean |pi - li| = {np.abs(realpi - li).mean():.2f}; mean |pi - R| = {np.abs(realpi - Rx).mean():.2f}; mean |pi - (R - sum li(x^rho))| = {np.abs(realpi - pred).mean():.2f}; corr(pi - R, zeros term) = {np.corrcoef(realpi - Rx, zsum)[0,1]:+.3f}, slope {np.polyfit(zsum, realpi - Rx, 1)[0]:.3f}")
for x in (1000, 10**4, 10**5, 10**6):
    i = np.searchsorted(xs, x) - 1
    print(f"    x={int(xs[i]):7d}: pi = {realpi[i]:6d}, li = {li[i]:9.1f}, R = {Rx[i]:9.1f}, R - sum li(x^rho) = {pred[i]:9.1f}")
print(f"[{time.time()-t0:.0f} s]")
