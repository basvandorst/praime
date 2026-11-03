import sys, time
import numpy as np
from math import log, sqrt
sys.path.insert(0, "/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins")
from explib import sieve, zeros, spectrum, main_term, hl_constant, C2

X = 10**8
L = 10**3
NZ = 300
NM = 100
gam = zeros(NZ + 1)
z = gam[:NZ]
mids = 0.5 * (gam[:NM] + gam[1:NM + 1])
ts = np.concatenate([z, mids])
main = main_term(X, L, ts)
u = log(X / L)
t0 = time.time()
print("E-D4b decoherence prediction for (p, p+h): ratio(h) = 1 + mean_j (1/ln(X/L)) int_L^X cos(gamma_j ln(1 + h/p)) dp/p")
lp = np.linspace(log(L), log(X), 4001)
p = np.exp(lp)
measured = {2: 1.999, 4: 1.830, 6: 1.793, 10: 1.805, 14: 1.925, 20: 1.842, 30: 1.773, 50: 1.559, 100: 1.541, 200: 1.858, 1000: 1.133}
print("     h   predicted   measured (300 zeros, +-0.17)")
for h in (2, 4, 6, 10, 14, 20, 30, 50, 100, 200, 1000):
    ph = np.log1p(h / p)
    c = np.array([np.trapezoid(np.cos(g * ph), lp) / u for g in z])
    print(f"  {h:5d}   {1 + c.mean():7.3f}     {measured[h]:.3f}")

isp = sieve(3 * X + 8)
primes = np.flatnonzero(isp[:X + 1])
primes = primes[primes >= L]
pf = primes.astype(np.float64)


def dip(vals, w):
    r = (spectrum(vals, w, ts) - main).real / u
    d = r[:NZ].mean() - r[NZ:].mean()
    e = sqrt((r[:NZ].std() / sqrt(NZ)) ** 2 + (r[NZ:].std() / sqrt(NM)) ** 2)
    return d, e, r


print("E-D3 scale pairs with the right singular series: (p, 2p+-1) constant 2C2; (p, 3p+-2) constant 2 * 2C2")
for name, mult, add, const in (("(p, 2p+1)", 2, 1, 2 * C2), ("(p, 2p-1)", 2, -1, 2 * C2), ("(p, 3p+2)", 3, 2, 4 * C2), ("(p, 3p-2)", 3, -2, 4 * C2)):
    v = pf[isp[mult * primes + add]]
    d, e, r = dip(v, np.log(v) * np.log(mult * v + add) / const)
    cz = np.cos(z * log(mult))
    a, b = np.polyfit(cz, -r[:NZ], 1)
    print(f"  {name}: n={len(v):7d}, weighted density {(np.log(v) * np.log(mult * v + add) / const).sum() / (X - L):.4f}, mean ratio {-d:5.3f} +- {e:.3f} (pred 1); fit {b:+.3f} {a:+.3f} cos(g ln {mult}) (pred 1 + {1/sqrt(mult):.3f} cos), corr {np.corrcoef(cz, -r[:NZ])[0,1]:+.3f}  [{time.time()-t0:.0f} s]", flush=True)
