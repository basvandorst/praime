import sys, time
import numpy as np
from math import log, sqrt, pi
from scipy.special import expi
sys.path.insert(0, "/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins")
from explib import sieve, primes_upto

X = 10**8
SEG = 10**7
t0 = time.time()
isp = sieve(X)
base = primes_upto(int(sqrt(X)) + 1)
mu = np.ones(X + 1, dtype=np.int8); mu[0] = 0
lam = np.ones(X + 1, dtype=np.int8); lam[0] = 0
for p in np.flatnonzero(isp):
    mu[p::p] *= -1
    q = p
    while q <= X:
        lam[q::q] *= -1
        q *= p
    if p <= X // p:
        mu[p * p::p * p] = 0
print(f"sieve, mu, lambda to {X:.0e}  [{time.time()-t0:.0f} s]", flush=True)
Lam = np.zeros(X + 1)
for p in base:
    q = p * p
    while q <= X:
        Lam[q] = log(p); q *= p
pr = np.flatnonzero(isp)
Lam[pr] = np.log(pr)
psi_off = theta_off = 0.0
pi_off = 0
M_off = L_off = 0
best = {"schoenfeld": (0, 0), "koch": (0, 0), "pili": (0, 0), "mertens": (0, 0), "polya": (0, 0)}
ratios = []
for lo in range(1, X + 1, SEG):
    hi = min(lo + SEG, X + 1)
    x = np.arange(lo, hi, dtype=np.float64)
    psi = psi_off + np.cumsum(Lam[lo:hi]); psi_off = psi[-1]
    th_inc = np.where(isp[lo:hi], np.log(np.maximum(x, 2)), 0.0)
    theta = theta_off + np.cumsum(th_inc); theta_off = theta[-1]
    pic = pi_off + np.cumsum(isp[lo:hi].astype(np.int64)); pi_off = pic[-1]
    M = M_off + np.cumsum(mu[lo:hi].astype(np.int64)); M_off = M[-1]
    Lv = L_off + np.cumsum(lam[lo:hi].astype(np.int64)); L_off = Lv[-1]
    m = x >= 1000
    xs, lx = x[m], np.log(x[m])
    r1 = np.abs(psi[m] - xs) / (np.sqrt(xs) * lx**2 / (8 * pi))
    r2 = np.abs(theta[m] - xs) / np.sqrt(xs)
    li = expi(lx) - expi(log(2))
    r3 = np.abs(pic[m] - li) / (np.sqrt(xs) * lx / (8 * pi))
    r4 = np.abs(M[m]) / np.sqrt(xs)
    r5 = Lv[m] / np.sqrt(xs)
    for key, r in (("schoenfeld", r1), ("koch", r2), ("pili", r3), ("mertens", r4), ("polya", r5)):
        i = np.argmax(r if key != "polya" else r)
        if r[i] > best[key][0]:
            best[key] = (float(r[i]), int(xs[i]))
    ratios.append((M[m][::997] / np.sqrt(xs[::997])))
    if lo == 1:
        print(f"  checks: pi(1e7) = {pic[10**7 - lo]} (664579), M(1e7) = {M[10**7 - lo]} (1037), L(1e7) = {Lv[10**7 - lo]}  [{time.time()-t0:.0f} s]", flush=True)
print("E-C3 Schoenfeld (RH <=> |psi(x) - x| < sqrt(x) ln^2 x / (8 pi) for x >= 73.2): max ratio for 1000 <= x <= 1e8: %.4f at x = %d" % best["schoenfeld"])
print("E-C3b von Koch form |theta(x) - x| / sqrt x: max %.4f at x = %d  (theta(x) < x throughout)" % best["koch"])
print("E-C4 |pi(x) - li(x)| < sqrt(x) ln x / (8 pi) (RH-equivalent, Schoenfeld): max ratio %.4f at x = %d" % best["pili"])
print("E-C5 Mertens |M(x)| / sqrt x: max %.4f at x = %d (known 0.5706 at 7.77e9 below 1e10; conjecture |M| < sqrt x is false far out)" % best["mertens"])
print("E-C5b Polya L(x)/sqrt x: max %.4f at x = %d (L(x) > 0 only in one stretch below 1e10, FACTS 2028)" % best["polya"])
allr = np.concatenate(ratios)
print(f"E-G1 M(x)/sqrt x sampled at {len(allr)} points in [1000, 1e8]: mean {allr.mean():+.4f}, sd {allr.std():.4f}, |max| {np.abs(allr).max():.4f}, kurtosis {((allr-allr.mean())**4).mean()/allr.var()**2:.2f}")
print("E-G2 Chowla: mean of lambda(n) lambda(n+h) and mu(n) mu(n+h) over n <= 1e8 (conjecture: -> 0 for lambda; mu(n)mu(n+h) -> 0 as well, while mu^2 pairs have density prod(1 - 2/p^2) = 0.3226):")
lf = lam[1:].astype(np.float64); mf = mu[1:].astype(np.float64)
for h in (1, 2, 3, 4, 6):
    cl = (lf[:-h] * lf[h:]).mean(); cm = (mf[:-h] * mf[h:]).mean(); sq = ((mf[:-h] != 0) & (mf[h:] != 0)).mean()
    print(f"  h={h}: <lambda lambda> {cl:+.5f}   <mu mu> {cm:+.5f}   P(both squarefree) {sq:.5f}   (1/sqrt N noise {1/sqrt(X):.5f})")
print(f"[{time.time()-t0:.0f} s]")
