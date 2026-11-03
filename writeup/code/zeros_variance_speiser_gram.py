import sys, time
import numpy as np
import mpmath as mp
from math import log, sqrt, pi, exp
sys.path.insert(0, "/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins")
from explib import zeros, lam_array, psi_from_zeros, GAMMA

t0 = time.time()
gam = zeros()
X = 2 * 10**6
Lam = lam_array(X)
psi = np.cumsum(Lam)
print("E-V1 variance of psi(x+H) - psi(x) - H over x in [1e6, 2e6] against Montgomery-Soundararajan H (ln(x/H) - gamma - ln 2pi) and against the same windows rebuilt from the zeros (x <= 2e6 is within the 1e5-zero resolution):")
xs = np.arange(10**6, 2 * 10**6 - 10**5, 977)
for H in (100, 1000, 10000, 100000):
    inc = psi[xs + H] - psi[xs] - H
    ms = H * (np.log(xs / H).mean() - GAMMA - log(2 * pi))
    print(f"  H={H:6d}: var(real) {inc.var():10.1f}   M-S {ms:10.1f}   ratio {inc.var()/ms:.3f}   mean {inc.mean():+.2f}", flush=True)
sub = xs[::7]
for H in (1000, 10000):
    pz = psi_from_zeros(np.concatenate([sub + 0.5, sub + H + 0.5]), gam)
    incz = pz[len(sub):] - pz[:len(sub)] - H
    incr = psi[sub + H] - psi[sub] - H
    print(f"  H={H:6d}: zeros-rebuilt windows: var {incz.var():10.1f}, corr with real {np.corrcoef(incz, incr)[0,1]:+.4f}, rms difference {sqrt(((incz-incr)**2).mean()):.2f}  [{time.time()-t0:.0f} s]", flush=True)

print("E-X2 Speiser: RH <=> zeta'(s) has no zeros with 0 < Re s < 1/2. Zeros of zeta' with 0 < Im s < 100 (found by Newton from a grid), and their real parts:")
mp.mp.dps = 15
found = []
for t in np.arange(2.0, 100.0, 0.5):
    for sig in (0.6, 1.0, 1.5):
        try:
            r = mp.findroot(lambda s: mp.zeta(s, derivative=1), mp.mpc(sig, t))
        except Exception:
            continue
        if 0 < r.imag < 100 and -0.5 < r.real < 3 and all(abs(r - f) > 1e-6 for f in found):
            found.append(r)
found.sort(key=lambda r: r.imag)
re = np.array([float(r.real) for r in found])
print(f"  {len(found)} zeros of zeta' found with Im < 100; min Re = {re.min():.4f}, max Re = {re.max():.4f}; count with Re < 1/2: {(re < 0.5).sum()}")
print("  first six:", ", ".join(f"{float(r.real):.4f}+{float(r.imag):.4f}i" for r in found[:6]))

print("E-X3 Gram's law: Gram points g_n solve theta(g_n) = n pi; the law says (-1)^n Z(g_n) > 0. Fraction of Gram intervals (g_n, g_{n+1}) holding exactly one zero, for the first 10000 zeros:")
def theta(t):
    return float(mp.siegeltheta(t))
gs = []
n = 0
t = 17.0
while n < 10000:
    t = float(mp.findroot(lambda u: mp.siegeltheta(u) - n * pi, t))
    gs.append(t); n += 1; t += 1.0
gs = np.array(gs)
counts = np.searchsorted(gam, gs[1:]) - np.searchsorted(gam, gs[:-1])
print(f"  exactly one zero: {(counts == 1).mean():.4f}; none: {(counts == 0).mean():.4f}; two: {(counts == 2).mean():.4f}; three: {(counts == 3).mean():.4f}  (literature: Gram's law holds for about 73 % of intervals, failures come in pairs)")
viol = np.flatnonzero(counts != 1)
print(f"  first violations at Gram indices: {viol[:8].tolist()} (first known failure near n = 126, t = 282.45)")
print(f"[{time.time()-t0:.0f} s]")
