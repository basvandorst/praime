import sys, time
import numpy as np
from math import log, sqrt
sys.path.insert(0, "/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins")
from explib import sieve, zeros, spectrum, main_term, hl_constant, C2

X = 10**8
NZ = 1000
NM = 300
gam = zeros(NZ + 1)
z = gam[:NZ]
mids = 0.5 * (gam[:NM] + gam[1:NM + 1])
ts = np.concatenate([z, mids])
t0 = time.time()
isp = sieve(2 * X + 8)
allp = np.flatnonzero(isp[:X + 1])
cst = {2: hl_constant([0, 2]), 4: hl_constant([0, 4]), 6: hl_constant([0, 6]), "t": hl_constant([0, 2, 6])}


def dip(vals, w, L):
    r = (spectrum(vals, w, ts) - main_term(X, L, ts)).real / log(X / L)
    d = r[:NZ].mean() - r[NZ:].mean()
    e = sqrt((r[:NZ].std() / sqrt(NZ)) ** 2 + (r[NZ:].std() / sqrt(NM)) ** 2)
    return d, e, r


for L in (10**3, 10**4):
    primes = allp[allp >= L]
    pf = primes.astype(np.float64)
    sub = pf[::7]
    dp, ep, _ = dip(sub, 7 * np.log(sub), L)
    print(f"L = {L}: prime dip (every 7th prime, weight 7 ln p) {dp:+.3f} +- {ep:.3f}  [{time.time()-t0:.0f} s]", flush=True)
    out = []
    for name, mask, w_fn, pred in (
        ("twins", isp[primes + 2], lambda v: np.log(v) ** 2 / cst[2], 2),
        ("cousins", isp[primes + 4], lambda v: np.log(v) ** 2 / cst[4], 2),
        ("sexy", isp[primes + 6], lambda v: np.log(v) ** 2 / cst[6], 2),
        ("triplets", isp[primes + 2] & isp[primes + 6], lambda v: np.log(v) ** 3 / cst["t"], 3),
        ("Sophie Germain", isp[2 * primes + 1], lambda v: np.log(v) * np.log(2 * v + 1) / cst[2], "1 + 0.707 cos"),
    ):
        v = pf[mask]
        d, e, r = dip(v, w_fn(v), L)
        line = f"  {name:15s} n={len(v):7d} ratio {d/dp:5.3f} +- {e/abs(dp):.3f}  (predicted {pred})"
        if name == "Sophie Germain":
            c = np.cos(z * log(2))
            a, b = np.polyfit(c, r[:NZ] / abs(dp), 1)
            line += f"; fit {b:+.3f} {a:+.3f} cos(g ln 2), corr {np.corrcoef(c, r[:NZ])[0,1]:+.3f} (predicted -1 -0.707 cos)"
        print(line + f"  [{time.time()-t0:.0f} s]", flush=True)
