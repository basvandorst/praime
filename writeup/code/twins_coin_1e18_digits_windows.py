import sys, time
import numpy as np
import gmpy2
from math import log, sqrt, exp
sys.path.insert(0, "/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins")
from explib import sieve, primes_upto, GAMMA, C2

t0 = time.time()
print("E-S1 twin coin at 1e18: survivors of the sieve to y = 1e18^0.4 = 1.6e7 in four windows of 2e6 at 1e18 + k*2.5e17, primality by BPSW; P(both prime)/(P(n) P(n+2)) and the single-coin rate")
y = int(1e18 ** 0.4)
ps = primes_upto(y)
rho_pred = 1 / (log(1e18) * np.prod(1 - 1 / ps[:2000].astype(np.float64)) * exp(sum(np.log1p(-1 / ps[2000:].astype(np.float64)))) * exp(GAMMA) * (1 + log(1 / 0.4 - 1)) / (1 / 0.4))
tot = {"n": 0, "a": 0, "b": 0, "c": 0}
for k in range(4):
    lo = 10**18 + k * 25 * 10**16
    W = 2 * 10**6
    comp = np.zeros(W + 2, dtype=bool)
    for p in ps:
        comp[(-lo) % p::p] = True
    S = ~comp
    T = np.flatnonzero(S[:W] & S[2:W + 2])
    a = np.array([gmpy2.is_prime(int(lo + t)) for t in T])
    b = np.array([gmpy2.is_prime(int(lo + t + 2)) for t in T])
    tot["n"] += len(T); tot["a"] += a.sum(); tot["b"] += b.sum(); tot["c"] += (a & b).sum()
    print(f"  window {k}: {len(T)} twin survivors, P(n prime) {a.mean():.4f}, P(n+2 prime) {b.mean():.4f}, both {(a&b).mean():.4f}  [{time.time()-t0:.0f} s]", flush=True)
n, a, b, c = tot["n"], tot["a"] / tot["n"], tot["b"] / tot["n"], tot["c"] / tot["n"]
print(f"  all: {n} twin survivors; rho measured {a:.4f} (model {rho_pred:.4f}); independence ratio c/(ab) = {c/(a*b):.4f} +- {sqrt((1-c)/(c*n)):.4f}")

print("E-S2 last digits of consecutive twin pairs (lesser members) to 1e8: transition matrix among last digits 1, 7, 9 against the product of marginals (the Lemke Oliver-Soundararajan bias for twins):")
isp = sieve(10**8 + 2)
pr = np.flatnonzero(isp[:10**8 + 1])
tw = pr[isp[pr + 2]]; tw = tw[tw > 5]
d = tw % 10
digits = (1, 7, 9)
M = np.zeros((3, 3))
for i, x in enumerate(digits):
    for j, y_ in enumerate(digits):
        M[i, j] = ((d[:-1] == x) & (d[1:] == y_)).sum()
M /= M.sum()
marg = M.sum(axis=1)
print("  P(next last digit | this last digit) / P(next last digit):")
for i, x in enumerate(digits):
    print(f"    from {x}: " + "  ".join(f"to {y_}: {M[i,j]/marg[i]/marg[j]:.3f}" for j, y_ in enumerate(digits)))
print(f"  diagonal share {np.trace(M):.4f} against {(marg**2).sum():.4f} for independence")

print("E-S3 twin counts in windows of length H near 1e8 (on [1e7, 1e8]): var/mean against Poisson (1) and the prime counts' values:")
tw_ind = np.zeros(10**8 + 1, dtype=bool); tw_ind[tw] = True
pr_ind = isp[:10**8 + 1]
for H in (10**3, 10**4, 10**5, 10**6):
    a_ = tw_ind[10**7:10**8].reshape(-1, H).sum(axis=1).astype(np.float64)
    b_ = pr_ind[10**7:10**8].reshape(-1, H).sum(axis=1).astype(np.float64)
    t = np.linspace(0, 1, len(a_))
    a_ = a_ - np.polyval(np.polyfit(t, a_, 1), t) + a_.mean()
    b_ = b_ - np.polyval(np.polyfit(t, b_, 1), t) + b_.mean()
    print(f"  H={H:7d}: twins var/mean {a_.var()/a_.mean():.3f} (mean {a_.mean():.1f});  primes var/mean {b_.var()/b_.mean():.3f}")
print(f"[{time.time()-t0:.0f} s]")
