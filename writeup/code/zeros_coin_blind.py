import sys, time
import numpy as np
from math import log, sqrt, exp
sys.path.insert(0, "/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins")
from explib import sieve, zeros, spectrum, main_term, primes_upto, GAMMA, C2

X = 10**7
L = 10**5
t0 = time.time()
isp = sieve(X + 2)
gam = zeros(201)
z = gam[:200]
mids = 0.5 * (gam[:60] + gam[1:61])
ts = np.concatenate([z, mids])
u = log(X / L)

print("E-S5 where do the zeros live, sieve layer or coin? Sieve to y = 100 over [1e5, 1e7]: survivors S, coin c(n) = 1_prime(n).")
y = 100
surv = np.ones(X + 3, dtype=bool)
for q in primes_upto(y):
    surv[::q] = False
idx = np.flatnonzero(surv[L:X + 1]) + L
nf = idx.astype(np.float64)
V = np.prod(1 - 1 / primes_upto(y).astype(np.float64))
rho = 1 / (np.log(nf) * V)
isprime = isp[idx]
print(f"  {len(idx)} survivors, prime share {isprime.mean():.4f}, mean rho {rho.mean():.4f}  [{time.time()-t0:.0f} s]")
main = main_term(X, L, ts)
S_prime = (spectrum(nf[isprime], np.log(nf[isprime]), ts) - main).real / u
S_sieve = (spectrum(nf, rho * np.log(nf), ts) - main).real / u
S_coin = (spectrum(nf, (isprime - rho) * np.log(nf), ts)).real / u
k = 200
print(f"  primes:            dips at zeros {S_prime[:k].mean():+.3f} +- {S_prime[:k].std()/sqrt(k):.3f}, midpoints {S_prime[k:].mean():+.3f}")
print(f"  sieve layer rho(n) over survivors: {S_sieve[:k].mean():+.3f} +- {S_sieve[:k].std()/sqrt(k):.3f}, midpoints {S_sieve[k:].mean():+.3f}   (carries the smooth term, no zeros expected)")
print(f"  coin (1_prime - rho) over survivors: {S_coin[:k].mean():+.3f} +- {S_coin[:k].std()/sqrt(k):.3f}, midpoints {S_coin[k:].mean():+.3f}   (all of the zeros expected here)")

print("E-S6 the same for twins: twin-survivor coin (1_twin - rho^2 C2tail) over twin survivors of the sieve to 100:")
tw_surv = surv[L:X + 1] & surv[L + 2:X + 3]
idx2 = np.flatnonzero(tw_surv) + L
nf2 = idx2.astype(np.float64)
V2 = 0.5 * np.prod(1 - 2 / primes_upto(y)[1:].astype(np.float64))
tail = np.prod([1 - 1 / (q - 1) ** 2 for q in primes_upto(10**6) if q > y])
rho2 = 2 * C2 / (np.log(nf2) ** 2 * V2)
istwin = isp[idx2] & isp[idx2 + 2]
print(f"  {len(idx2)} twin survivors, twin share {istwin.mean():.5f}, mean rho2 {rho2.mean():.5f}")
T_coin = (spectrum(nf2, (istwin - rho2) * np.log(nf2) ** 2 / (2 * C2), ts)).real / u
print(f"  twin coin spectrum: dips at zeros {T_coin[:k].mean():+.3f} +- {T_coin[:k].std()/sqrt(k):.3f}, midpoints {T_coin[k:].mean():+.3f}   (expected about twice the prime coin)  [{time.time()-t0:.0f} s]")

print("E-A8 blind zero finding from the primes alone (X = 1e7, p >= 1e3): scan t in [10, 60] with step 0.02, refine minima by a parabola, compare with Odlyzko:")
LB = 10**3
pr = np.flatnonzero(isp[:X + 1]); pr = pr[pr >= LB].astype(np.float64)
tg = np.arange(10.0, 60.0, 0.02)
ub = log(X / LB)
Sg = (spectrum(pr, np.log(pr), tg) - main_term(X, LB, tg)).real / ub
found = []
for i in range(1, len(tg) - 1):
    if Sg[i] < Sg[i - 1] and Sg[i] < Sg[i + 1] and Sg[i] < -0.5:
        a, b, c = Sg[i - 1], Sg[i], Sg[i + 1]
        off = 0.5 * (a - c) / (a - 2 * b + c)
        found.append(tg[i] + off * 0.02)
found = np.array(found)
true = gam[(gam > 10) & (gam < 60)]
print(f"  {len(found)} dips found, {len(true)} true zeros in (10, 60)")
match = np.array([np.abs(found - g).min() for g in true])
print(f"  distance of each true zero to the nearest dip: max {match.max():.4f}, rms {sqrt((match**2).mean()):.4f}; spurious dips: {sum(np.abs(true - f).min() > 0.1 for f in found)}")
print("  first eight recovered:", ", ".join(f"{f:.4f}" for f in found[:8]))
print("  Odlyzko:               ", ", ".join(f"{g:.4f}" for g in true[:8]))
print(f"[{time.time()-t0:.0f} s]")
