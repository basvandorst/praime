import sys, time
import numpy as np
from math import log, sqrt
sys.path.insert(0, "/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins")
from explib import sieve, zeros, spectrum, main_term, C2

X = 10**8
L = 10**3
NZ = 1000
NM = 300
PRIME_DIP = -1.008
gam = zeros(NZ + 1)
z = gam[:NZ]
mids = 0.5 * (gam[:NM] + gam[1:NM + 1])
ts = np.concatenate([z, mids])
main = main_term(X, L, ts)
u = log(X / L)
t0 = time.time()
isp = sieve(X + 2)
pi = np.cumsum(isp[:X + 1], dtype=np.int32)
primes = np.flatnonzero(isp[:X + 1])
primes = primes[primes >= L]
pf = primes.astype(np.float64)
lnp = np.log(pf)
w = (pf * 1e-3).astype(np.int64)
local = (pi[np.minimum(primes + w, X)] - pi[np.maximum(primes - w, 0)]).astype(np.float64) / (2 * w / lnp)
print(f"local density D(p) over windows +-p/1000: mean {local.mean():.4f}, sd {local.std():.4f}  [{time.time()-t0:.0f} s]")
rng = np.random.default_rng(11)


def ratio(vals, wts):
    r = (spectrum(vals, wts, ts) - main).real / u
    d = r[:NZ].mean() - r[NZ:].mean()
    return d / PRIME_DIP


real = primes[isp[primes + 2]].astype(np.float64)
print(f"real twins ({len(real)}, p >= {L}): ratio {ratio(real, np.log(real) ** 2 / (2 * C2)):.3f}  [{time.time()-t0:.0f} s]")
r1s, r2s = [], []
for trial in range(4):
    k1 = rng.random(len(pf)) < 2 * C2 / lnp
    r1 = ratio(pf[k1], lnp[k1] ** 2 / (2 * C2))
    prob = np.minimum(1.0, 2 * C2 / lnp * local)
    k2 = rng.random(len(pf)) < prob
    r2 = ratio(pf[k2], lnp[k2] ** 2 / (2 * C2))
    r1s.append(r1); r2s.append(r2)
    print(f"trial {trial}: coin partner ({k1.sum()}) ratio {r1:.3f} (truth 1);  local-density partner ({k2.sum()}) ratio {r2:.3f} (truth 2)  [{time.time()-t0:.0f} s]", flush=True)
r1s, r2s = np.array(r1s), np.array(r2s)
print(f"means: coin {r1s.mean():.3f} +- {r1s.std(ddof=1)/2:.3f} (truth 1); local-density {r2s.mean():.3f} +- {r2s.std(ddof=1)/2:.3f} (truth 2)")
