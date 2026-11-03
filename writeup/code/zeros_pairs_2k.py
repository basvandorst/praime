import sys, time
import numpy as np
from math import log, sqrt
sys.path.insert(0, "/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins")
from explib import sieve, zeros, spectrum, main_term, hl_constant, C2

X = 10**8
L = 10**3
NZ = 300
NM = 100
PRIME_DIP = -1.0
gam = zeros(NZ + 1)
z = gam[:NZ]
mids = 0.5 * (gam[:NM] + gam[1:NM + 1])
ts = np.concatenate([z, mids])
main = main_term(X, L, ts)
u = log(X / L)
t0 = time.time()
isp = sieve(3 * X + 8)
primes = np.flatnonzero(isp[:X + 1])
primes = primes[primes >= L]
pf = primes.astype(np.float64)


def dip(vals, w):
    r = (spectrum(vals, w, ts) - main).real / u
    d = r[:NZ].mean() - r[NZ:].mean()
    e = sqrt((r[:NZ].std() / sqrt(NZ)) ** 2 + (r[NZ:].std() / sqrt(NM)) ** 2)
    return d, e, r


print("E-D4 pairs (p, p + 2k) at the first 300 zeros, X = 1e8, p >= 1000; circle method: ratio 2 for every k")
for k in (1, 2, 3, 5, 7, 10, 15, 25, 50, 100, 500):
    h = 2 * k
    v = pf[isp[primes + h]]
    c = hl_constant([0, h])
    d, e, _ = dip(v, np.log(v) * np.log(v + h) / c)
    print(f"  h={h:4d}: n={len(v):7d}, S(h)={c:.4f}, ratio {-d:5.3f} +- {e:.3f}  [{time.time()-t0:.0f} s]", flush=True)
print("E-D3 pairs at different scales: (p, 2p+1) predicts 1 + 0.707 cos(g ln 2); (p, 2p-1) the same; (p, 3p+2) predicts 1 + 0.577 cos(g ln 3)")
for name, mult, add, base in (("(p, 2p+1)", 2, 1, 2), ("(p, 2p-1)", 2, -1, 2), ("(p, 3p+2)", 3, 2, 3), ("(p, 3p-2)", 3, -2, 3)):
    v = pf[isp[mult * primes + add]]
    c = hl_constant([0, 1]) if base == 2 else hl_constant([0, 2])
    w = np.log(v) * np.log(mult * v + add) / c
    d, e, r = dip(v, w)
    cz = np.cos(z * log(base))
    a, b = np.polyfit(cz, -r[:NZ], 1)
    print(f"  {name}: n={len(v):7d}, mean ratio {-d:5.3f} +- {e:.3f} (pred 1); fit {b:+.3f} {a:+.3f} cos(g ln {base}) (pred 1 + {1/sqrt(base):.3f} cos), corr {np.corrcoef(cz, -r[:NZ])[0,1]:+.3f}  [{time.time()-t0:.0f} s]", flush=True)
