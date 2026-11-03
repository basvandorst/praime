"""Odlyzko's zeros: curl -sL -o zeros1 http://www.dtc.umn.edu/~odlyzko/zeta_tables/zeros1 (first 100000 zeros, 1.8 MB), kept outside the repo.
Usage: python3 -I zeros_twin_spectrum.py <path to zeros1>
"""
import sys, time
import numpy as np
from math import log, sqrt

X = 10**8
C2 = 0.6601618158468696
NZ = 300
NM = 100
gam = np.loadtxt(sys.argv[1], max_rows=NZ + 1)
z = gam[:NZ]
mids = 0.5 * (gam[:NM] + gam[1:NM + 1])
ts = np.concatenate([z, mids])
edges = [2, 10**2, 10**4, 10**6, 10**8]


def small_primes(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if s[i]: s[i * i::i] = False
    return s


def spectrum(values, weights, ts, chunk=20000):
    acc = np.zeros(len(ts), dtype=np.complex128); lv = np.log(values); w = weights / np.sqrt(values)
    for i in range(0, len(values), chunk):
        acc += (w[None, i:i + chunk] * np.exp(-1j * np.outer(ts, lv[i:i + chunk]))).sum(axis=1)
    return acc


t0 = time.time()
isp = small_primes(X + 2)
primes = np.flatnonzero(isp[:X + 1]); twins = primes[isp[primes + 2]]
pr = primes.astype(np.float64); tw = twins.astype(np.float64)
SP, ST = [], []
for a, b in zip(edges[:-1], edges[1:]):
    m = (pr >= a) & (pr < b); SP.append(spectrum(pr[m], np.log(pr[m]), ts))
    m = (tw >= a) & (tw < b); ST.append(spectrum(tw[m], np.log(tw[m]) ** 2 / (2 * C2), ts))
    print(f"  band [{a:.0e}, {b:.0e}) done [{time.time()-t0:.0f} s]", flush=True)
print("dips in units of ln(X/L), primes without prime powers; mean over the first 300 zeros, baseline over 100 midpoints")
print("  lower cut L   primes: zeros   midpoints     twins: zeros    midpoints     ratio of (zero - midpoint)")
for k, L in enumerate(edges[:-1]):
    main = (X ** (0.5 - 1j * ts) - L ** (0.5 - 1j * ts)) / (0.5 - 1j * ts)
    u = log(X / L)
    rp = (sum(SP[k:]) - main).real / u
    rt = (sum(ST[k:]) - main).real / u
    dp = rp[:NZ].mean() - rp[NZ:].mean(); dt = rt[:NZ].mean() - rt[NZ:].mean()
    err = sqrt((rt[:NZ].std() / sqrt(NZ)) ** 2 + (rt[NZ:].std() / sqrt(NM)) ** 2) / abs(dp)
    print(f"  {L:8.0e}      {rp[:NZ].mean():+.3f} +- {rp[:NZ].std()/sqrt(NZ):.3f}  {rp[NZ:].mean():+.3f}        {rt[:NZ].mean():+.3f} +- {rt[:NZ].std()/sqrt(NZ):.3f}  {rt[NZ:].mean():+.3f}        {dt/dp:.2f} +- {err:.2f}")
print("band-only dips (each band normalised by its own ln(b/a)):")
for k, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
    main = (b ** (0.5 - 1j * ts) - a ** (0.5 - 1j * ts)) / (0.5 - 1j * ts)
    u = log(b / a)
    rp = (SP[k] - main).real / u; rt = (ST[k] - main).real / u
    print(f"  [{a:.0e},{b:.0e})  primes {rp[:NZ].mean():+.3f} (mid {rp[NZ:].mean():+.3f})   twins {rt[:NZ].mean():+.3f} +- {rt[:NZ].std()/sqrt(NZ):.3f} (mid {rt[NZ:].mean():+.3f})   ratio {(rt[:NZ].mean()-rt[NZ:].mean())/(rp[:NZ].mean()-rp[NZ:].mean()):.2f}")
