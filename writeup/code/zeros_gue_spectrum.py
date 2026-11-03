"""Odlyzko's zeros: curl -sL -o zeros1 http://www.dtc.umn.edu/~odlyzko/zeta_tables/zeros1 (first 100000 zeros, 1.8 MB), kept outside the repo.
Usage: python3 -I zeros_gue_spectrum.py <path to zeros1>
"""
import sys, time
import numpy as np
from math import log, sqrt, pi

ZEROS = sys.argv[1]
X = 10**8
NZ = 100
C2 = 0.6601618158468696
gam = np.loadtxt(ZEROS, max_rows=100000)


def small_primes(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return s


t0 = time.time()
isp = small_primes(X + 2)
primes = np.flatnonzero(isp[:X + 1])
twins = primes[isp[primes + 2]]
print(f"primes <= {X:.0e}: {len(primes)}, twin pairs: {len(twins)}  [{time.time()-t0:.0f} s]")

ns, lam = [primes.astype(np.float64)], [np.log(primes.astype(np.float64))]
for p in primes[primes <= int(sqrt(X))]:
    q = p * p
    while q <= X:
        ns.append(np.array([q], dtype=np.float64)); lam.append(np.array([log(p)]))
        q *= p
n_all = np.concatenate(ns); lam_all = np.concatenate(lam)

z = gam[:NZ]
mids = 0.5 * (gam[:NZ] + gam[1:NZ + 1])
ts = np.concatenate([z, mids])


def spectrum(values, weights, ts, chunk=20000):
    acc = np.zeros(len(ts), dtype=np.complex128)
    lv = np.log(values)
    w = weights / np.sqrt(values)
    for i in range(0, len(values), chunk):
        l = lv[i:i + chunk]; ww = w[i:i + chunk]
        acc += (ww[None, :] * np.exp(-1j * np.outer(ts, l))).sum(axis=1)
    return acc


main = X ** (0.5 - 1j * ts) / (0.5 - 1j * ts)
SP = spectrum(n_all, lam_all, ts) - main
tw = twins.astype(np.float64)
ST = spectrum(tw, np.log(tw) ** 2 / (2 * C2), ts) - main
lnX = log(X)
print(f"spectra computed [{time.time()-t0:.0f} s]; depth unit ln X = {lnX:.2f}")
rp, rt = SP.real / lnX, ST.real / lnX
print()
print("Re S / ln X at the first 100 zeros and at the 100 midpoints between them (mean +- standard error):")
print(f"  primes, at zeros   {rp[:NZ].mean():+.3f} +- {rp[:NZ].std()/sqrt(NZ):.3f}      at midpoints {rp[NZ:].mean():+.3f} +- {rp[NZ:].std()/sqrt(NZ):.3f}")
print(f"  twins,  at zeros   {rt[:NZ].mean():+.3f} +- {rt[:NZ].std()/sqrt(NZ):.3f}      at midpoints {rt[NZ:].mean():+.3f} +- {rt[NZ:].std()/sqrt(NZ):.3f}")
print(f"  twins/primes depth ratio at zeros: {rt[:NZ].mean()/rp[:NZ].mean():.2f}  (1 = one coin inherits the density, 2 = both members do)")
print("  first ten zeros:  gamma      primes   twins")
for j in range(10):
    print(f"                 {z[j]:9.4f}   {rp[j]:+.3f}   {rt[j]:+.3f}")

print()
print("GUE statistics of the zeros 20000..100000 (unfolded spacings s = dgamma ln(gamma/2pi)/2pi):")
g = gam[20000:]
s = np.diff(g) * np.log(g[:-1] / (2 * pi)) / (2 * pi)
print(f"  mean s = {s.mean():.4f} (1), var s = {s.var():.4f} (GUE 0.178, Poisson 1.000)")
edges = np.array([0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 4.0])
hist = np.histogram(s, bins=edges)[0] / len(s)
from scipy.integrate import quad
wig = lambda x: 32 / pi**2 * x * x * np.exp(-4 * x * x / pi)
gue = np.array([quad(wig, a, b)[0] for a, b in zip(edges[:-1], edges[1:])])
poi = np.exp(-edges[:-1]) - np.exp(-edges[1:])
print("  bin          zeros    GUE surmise   Poisson")
for a, b, h, w_, p_ in zip(edges[:-1], edges[1:], hist, gue, poi):
    print(f"  [{a:4.2f},{b:4.2f})  {h:.4f}   {w_:.4f}        {p_:.4f}")
