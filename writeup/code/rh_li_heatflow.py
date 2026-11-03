import sys, time
import numpy as np
from math import log, pi, sqrt

GAMMA = 0.5772156649015329
gam = np.loadtxt(sys.argv[1], max_rows=100000)
rho = 0.5 + 1j * gam
t0 = time.time()
print("Li's criterion: RH <=> lambda_n = sum_rho [1 - (1 - 1/rho)^n] >= 0 for all n >= 1")
ns = np.array([1, 2, 3, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000])
one_minus = 1 - 1 / rho
lam = []
for n in ns:
    v = 2 * (1 - one_minus ** n).real.sum()
    lam.append(v)
lam = np.array(lam)
T = gam[-1]
tail = ns.astype(float) ** 2 * (log(T / (2 * pi)) + 1) / (2 * pi * T)
lam_c = lam + tail
smooth = ns / 2 * (np.log(ns) - 1 - log(2 * pi) + GAMMA)
print("     n      lambda_n (1e5 zeros + tail)   tail added   smooth (n/2)(ln n - 1 - ln 2pi + gamma)   difference   diff/(sqrt n ln n)")
for n, l, tl, sm in zip(ns, lam_c, tail, smooth):
    print(f"  {n:6d}   {l:16.4f}   {tl:10.3f}   {sm:16.4f}   {l-sm:+10.4f}   {(l-sm)/(sqrt(n)*log(n)) if n > 1 else float('nan'):+8.3f}")
print("  under RH the difference is O(sqrt n ln n) (Lagarias 2007); a zero off the line would make it grow like a power of n and lambda_n eventually negative")
print(f"  all positive: {bool((lam_c > 0).all())}; lambda_1 exact = 1 + gamma/2 - ln(4 pi)/2 = {1 + GAMMA/2 - log(4*pi)/2:.6f}  [{time.time()-t0:.0f} s]")

print()
print("de Bruijn-Newman: RH <=> Lambda <= 0, and Lambda >= 0 is proved (Rodgers-Tao 2018). Under the heat flow the zeros move by dx_j/dt = sum_k 2/(x_j - x_k);")
print("running it backwards from t = 0, the closest pairs collide first; the collision time of the closest pair in the first 100000 zeros is a (weak) lower bound for Lambda.")
d = np.diff(gam)
dens = np.log(gam[:-1] / (2 * pi)) / (2 * pi)
norm = d * dens
order = np.argsort(norm)[:5]
print("  closest pairs (normalised spacing s = gap * density):")
for i in order:
    print(f"    gamma = {gam[i]:.6f}, {gam[i+1]:.6f}: gap {d[i]:.6f}, s = {norm[i]:.4f}")
i = order[0]
x = gam.copy()
t, dt = 0.0, -1e-7
j = i
while True:
    lo, hi = max(0, j - 600), min(len(x), j + 602)
    xs = x[lo:hi]
    diff = xs[:, None] - xs[None, :]
    np.fill_diagonal(diff, np.inf)
    v = (2 / diff).sum(axis=1)
    gap = x[j + 1] - x[j]
    step = min(abs(dt), 0.05 * gap / max(abs(v[j + 1 - lo] - v[j - lo]), 1e-9))
    if gap < 1e-5 or t < -1e-2:
        break
    x[lo:hi] -= step * v
    t -= step
print(f"  backward flow of the closest pair (gamma {gam[i]:.4f}, {gam[i+1]:.4f}) with 1200 neighbours: collision near t = {t:.3e}  (isolated-pair estimate -gap^2/8 = {-d[i]**2/8:.3e})")
print(f"  so these zeros alone give Lambda > {t:.1e}; the literature's Lehmer-pair bound Lambda > -1.15e-11 (Saouter-Gourdon-Demichel 2011) needs pairs at height 1e21, and Rodgers-Tao give Lambda >= 0 exactly.")
