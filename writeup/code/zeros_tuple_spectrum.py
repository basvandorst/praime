import sys, time
import numpy as np
from math import log, sqrt

X = 10**8
L = 10**2
NZ = 300
NM = 100
C2 = 0.6601618158468696
gam = np.loadtxt(sys.argv[1], max_rows=NZ + 1)
z = gam[:NZ]
mids = 0.5 * (gam[:NM] + gam[1:NM + 1])
ts = np.concatenate([z, mids])


def small_primes(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if s[i]: s[i * i::i] = False
    return s


def hl_constant(offsets, qmax=10**6):
    k = len(offsets)
    ps = np.flatnonzero(small_primes(qmax))
    c = 1.0
    for q in ps:
        nu = len(set(o % q for o in offsets))
        c *= (1 - nu / q) / (1 - 1 / q) ** k
    return c


def spectrum(values, weights, ts, chunk=20000):
    acc = np.zeros(len(ts), dtype=np.complex128); lv = np.log(values); w = weights / np.sqrt(values)
    for i in range(0, len(values), chunk):
        acc += (w[None, i:i + chunk] * np.exp(-1j * np.outer(ts, lv[i:i + chunk]))).sum(axis=1)
    return acc


t0 = time.time()
isp = small_primes(2 * X + 8)
primes = np.flatnonzero(isp[:X + 1])
primes = primes[primes >= L].astype(np.float64)
pi = primes.astype(np.int64)
sets = {
    "primes (ln p)": (primes, np.log(primes), 1.0),
    "twins (p, p+2)": (primes[isp[pi + 2]], None, hl_constant([0, 2])),
    "cousins (p, p+4)": (primes[isp[pi + 4]], None, hl_constant([0, 4])),
    "sexy (p, p+6)": (primes[isp[pi + 6]], None, hl_constant([0, 6])),
    "triplets (p, p+2, p+6)": (primes[isp[pi + 2] & isp[pi + 6]], None, hl_constant([0, 2, 6])),
    "Sophie Germain (p, 2p+1)": (primes[isp[2 * pi + 1]], None, None),
}
main = (X ** (0.5 - 1j * ts) - L ** (0.5 - 1j * ts)) / (0.5 - 1j * ts)
u = log(X / L)
print(f"X = {X:.0e}, primes from {L}, {len(primes)} primes; dips in units of ln(X/L) = {u:.2f}; first {NZ} zeros, {NM} midpoints")
print(f"HL constants: twins {sets['twins (p, p+2)'][2]:.5f}, cousins {sets['cousins (p, p+4)'][2]:.5f}, sexy {sets['sexy (p, p+6)'][2]:.5f}, triplets {sets['triplets (p, p+2, p+6)'][2]:.5f}")
res = {}
for name, (vals, w, const) in sets.items():
    if name.startswith("primes"):
        pass
    elif name.startswith("Sophie"):
        const = hl_constant([0, 2])
        w = np.log(vals) * np.log(2 * vals + 1) / const
    elif name.startswith("triplets"):
        w = np.log(vals) ** 3 / const
    else:
        w = np.log(vals) ** 2 / const
    r = (spectrum(vals, w, ts) - main).real / u
    res[name] = r
    print(f"  {name:26s} n={len(vals):8d}  at zeros {r[:NZ].mean():+.3f} +- {r[:NZ].std()/sqrt(NZ):.3f}   midpoints {r[NZ:].mean():+.3f} +- {r[NZ:].std()/sqrt(NM):.3f}  [{time.time()-t0:.0f} s]", flush=True)
rp = res["primes (ln p)"]
dp = rp[:NZ].mean() - rp[NZ:].mean()
print()
print("depth ratios to the prime dip (zero minus midpoint), with the circle-method prediction:")
for name, pred in [("twins (p, p+2)", 2), ("cousins (p, p+4)", 2), ("sexy (p, p+6)", 2), ("triplets (p, p+2, p+6)", 3), ("Sophie Germain (p, 2p+1)", "1 + 0.707 cos(g ln 2), mean 1")]:
    r = res[name]
    d = r[:NZ].mean() - r[NZ:].mean()
    err = sqrt((r[:NZ].std() / sqrt(NZ)) ** 2 + (r[NZ:].std() / sqrt(NM)) ** 2) / abs(dp)
    print(f"  {name:26s} {d/dp:5.2f} +- {err:.2f}   predicted {pred}")
sg = res["Sophie Germain (p, 2p+1)"][:NZ]
c = np.cos(z * log(2))
a, b = np.polyfit(c, sg, 1)
print()
print(f"Sophie Germain per-zero pattern: fit dip = {b:+.3f} + {a:+.3f} cos(gamma ln 2); predicted -1.000 - 0.707 cos; corr(dip, cos) = {np.corrcoef(c, sg)[0,1]:+.3f}")
tw = res["twins (p, p+2)"][:NZ]
print(f"control: twins against cos(gamma ln 2): slope {np.polyfit(c, tw, 1)[0]:+.3f}, corr {np.corrcoef(c, tw)[0,1]:+.3f}")
