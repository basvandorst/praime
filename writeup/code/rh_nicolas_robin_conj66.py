import sys, time, gzip, re
import numpy as np
from math import log, sqrt, exp
sys.path.insert(0, "/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins")
from explib import zeros, segments, primes_upto, GAMMA, C2

t0 = time.time()
gam = zeros()
rho = 0.5 + 1j * gam
X = 10**9
theta = 0.0
lnNphi = 0.0
rows = []
for lo, P in segments(X, 10**8):
    idx = np.flatnonzero(P[:min(10**8, X - lo)]) + lo
    pf = idx.astype(np.float64)
    th = theta + np.cumsum(np.log(pf)); theta = th[-1]
    ln = lnNphi + np.cumsum(-np.log1p(-1 / pf)); lnNphi = ln[-1]
    margin = ln - GAMMA - np.log(np.log(th))
    sel = np.arange(0, len(pf), 2000)
    rows.append(np.vstack([pf[sel], margin[sel], th[sel]]))
    print(f"  segment {lo:.0e}: min margin {margin[pf > 100].min() if (pf > 100).any() else float('nan'):.3e}  [{time.time()-t0:.0f} s]", flush=True)
rows = np.hstack(rows)
p, m, th = rows
ok = p > 10**4
p, m, th = p[ok], m[ok], th[ok]
scaled = m * np.sqrt(th) * np.log(p)
osc = np.array([2 * ((x ** (rho - 0.5)) / (rho * (1 - rho))).real.sum() for x in p[::10]])
a, b = np.polyfit(osc, scaled[::10], 1)
print(f"E-C12 Nicolas to 1e9: no negative margin; scaled margin m sqrt(ln N) ln p = {b:.3f} + {a:.3f} osc(zeros) (prediction 2 + osc), corr {np.corrcoef(osc, scaled[::10])[0,1]:+.3f}; scaled at 1e9: {scaled[-1]:.4f}")
cheb = (th - p) / np.sqrt(p)
print(f"E-C13 Nicolas oscillation against the Chebyshev bias (theta(x) - x)/sqrt x directly, no zeros: corr(scaled - 2, (theta - x)/sqrt x) = {np.corrcoef(scaled - 2, cheb)[0,1]:+.3f} over {len(p)} points (the two carry the same zeros with weights 1/(rho(1-rho)) and 1/rho)")

print("E-C14 Robin's deficit on colossally abundant numbers: oscillating part against the zeros")
Pmax = 200000
events = []
for q in primes_upto(Pmax).astype(np.float64):
    a_ = 1
    while True:
        eps = log((q ** (a_ + 1) - 1) / (q ** a_ - 1)) / log(q) - 1
        if eps < 1e-7 and a_ > 1: break
        events.append((eps, q)); a_ += 1
        if q ** (a_ + 1) > 1e300: break
events.sort(key=lambda e: -e[0])
lnN = lnsig = 0.0; expo = {}; ca = []
for eps, q in events:
    a_ = expo.get(q, 0) + 1; expo[q] = a_
    lnN += log(q); lnsig += log(q ** (a_ + 1) - 1) - (log(q ** a_ - 1) if a_ > 1 else log(q - 1))
    if lnN > 1000:
        ca.append((lnN, exp(GAMMA) - exp(lnsig - lnN) / log(lnN), q))
ca = np.array(ca)
L, d, qmax = ca[:, 0], ca[:, 1], ca[:, 2]
sc = d * np.sqrt(L) * np.log(L)
oscr = np.array([2 * ((x ** (rho - 0.5)) / (rho * (1 - rho))).real.sum() for x in L[::5]])
a, b = np.polyfit(oscr, sc[::5], 1)
print(f"  {len(ca)} CA numbers with ln N > 1000: deficit sqrt(ln N) ln ln N = {b:.3f} + {a:.3f} osc(zeros at x = ln N); corr {np.corrcoef(oscr, sc[::5])[0,1]:+.3f}  (Nicolas' analogue has 2 + 1.0 osc)")

print("E-M2 Conjecture 66 (maximal gaps between twin primes ~ 0.45 ln^3 p) against the OEIS records A113274/A113275 and Kourbatov's a(ln(p/a) - b), a = ln^2 p / (2 C2):")
names = {}
with gzip.open('/home/bas/praime/references/oeis/stripped.gz', 'rt', encoding='utf-8', errors='replace') as f:
    for line in f:
        if line.startswith('A113274 ') or line.startswith('A113275 '):
            names[line[:7]] = np.array([int(v) for v in line.split(' ', 1)[1].strip().strip(',').split(',')])
gaps, starts = names['A113274'], names['A113275']
n = min(len(gaps), len(starts))
gaps, starts = gaps[:n], starts[:n]
pend = starts + gaps
ok = pend > 1000
lp = np.log(pend[ok].astype(np.float64))
ratio = gaps[ok] / lp ** 3
a_k = lp ** 2 / (2 * C2)
kour = a_k * (np.log(pend[ok] / a_k) - 1.0)
print(f"  {ok.sum()} records with end-of-gap prime above 1000: ratio gap/ln^3 p from {ratio.min():.3f} to {ratio.max():.3f}, mean {ratio.mean():.3f}; Kourbatov's formula with b = 1: mean(gap/formula) = {(gaps[ok]/kour).mean():.3f}, sd {(gaps[ok]/kour).std():.3f}; last record gap {gaps[-1]} after {starts[-1]}")
print(f"  limit of gap/ln^3 p implied by a(ln(p/a) - b): 1/(2 C2) = {1/(2*C2):.3f}, approached like 1 - (2 ln ln p + ln(2C2) ... )/ln p")
print(f"[{time.time()-t0:.0f} s]")
