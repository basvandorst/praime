import sys, time
import numpy as np
from math import log, exp, sqrt

GAMMA = 0.5772156649015329
P = int(float(sys.argv[1])) if len(sys.argv) > 1 else 300000
t0 = time.time()
s = np.ones(P + 1, dtype=bool); s[:2] = False
for i in range(2, int(P**0.5) + 1):
    if s[i]: s[i * i::i] = False
primes = np.flatnonzero(s).astype(np.float64)
events = []
for p in primes:
    a = 1
    while True:
        eps = log((p ** (a + 1) - 1) / (p ** a - 1)) / log(p) - 1
        if eps < 1e-7 and a > 1:
            break
        events.append((eps, p))
        a += 1
        if p ** (a + 1) > 1e300:
            break
events.sort(key=lambda e: -e[0])
print(f"{len(events)} colossally abundant steps from primes to {P}  [{time.time()-t0:.0f} s]")
lnN = 0.0
lnsig = 0.0
expo = {}
rows = []
for k, (eps, p) in enumerate(events):
    a = expo.get(p, 0) + 1
    expo[p] = a
    lnN += log(p)
    lnsig += log(p ** (a + 1) - 1) - (log(p ** a - 1) if a > 1 else log(p - 1))
    if lnN < 20:
        continue
    ratio = exp(lnsig - lnN) / log(lnN)
    deficit = exp(GAMMA) - ratio
    rows.append((lnN, deficit, p))
rows = np.array(rows)
lnN, deficit, pmax = rows[:, 0], rows[:, 1], rows[:, 2]
print("Robin (1984): RH <=> sigma(n) < e^gamma n ln ln n for all n > 5040; the colossally abundant numbers are where sigma(n)/(n ln ln n) is largest.")
print(f"  CA numbers with ln N in [20, {lnN[-1]:.3e}]: smallest deficit e^gamma - sigma/(N ln ln N) = {deficit.min():.3e} at ln N = {lnN[deficit.argmin()]:.1f}; negative values: {(deficit < 0).sum()}")
print("     ln N        largest p    deficit      deficit*sqrt(ln N)   deficit*sqrt(ln N)*ln ln N   deficit*sqrt(ln N)*(ln ln N)^2")
for target in (1e2, 1e3, 1e4, 1e5, 3e5):
    i = np.searchsorted(lnN, target)
    if i >= len(lnN): break
    L = lnN[i]; d = deficit[i]
    print(f"  {L:10.1f}  {int(pmax[i]):10d}   {d:.4e}   {d*sqrt(L):12.4f}   {d*sqrt(L)*log(L):14.4f}   {d*sqrt(L)*log(L)**2:16.4f}")
m = lnN > 1e3
a, b = np.polyfit(np.log(lnN[m]), np.log(deficit[m]), 1)
print(f"  power-law fit of the deficit for ln N > 1e3: deficit ~ (ln N)^({a:.3f});  RH-type scaling would be -1/2 times slowly varying factors")
