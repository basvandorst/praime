import sys, time
import numpy as np
import mpmath as mp
from math import log, sqrt, pi
sys.path.insert(0, "/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins")
from explib import zeros, psi_from_zeros, lam_array, sieve, primes_upto

gam = zeros()
t0 = time.time()
print("E-A1 truncation: max |Lambda(n) - increment of psi from N zeros| for n <= 2000, 10000 against N")
N = 10000
Lam = lam_array(N)[2:]
xs = np.arange(1.5, N + 1.0, 1.0)
for nz in (1000, 3000, 10000, 30000, 100000):
    psi = psi_from_zeros(xs, gam[:nz])
    err = np.diff(psi) - Lam
    T = gam[nz - 1]
    print(f"  N={nz:6d} (T={T:9.1f}): max|err| n<=2000 {np.abs(err[:1999]).max():.3f}, n<=10000 {np.abs(err).max():.3f};  n/T at the max: {np.argmax(np.abs(err))/T:.3f}  [{time.time()-t0:.0f} s]", flush=True)

print("E-A5 Mertens M(x) from zeros: M(x) ~ sum_rho x^rho/(rho zeta'(rho)) - 2 (first 1000 zeros, zeta' by mpmath)")
mp.mp.dps = 15
nz = 1000
dz = np.array([complex(mp.zeta(mp.mpc(0.5, g), derivative=1)) for g in gam[:nz]])
rho = 0.5 + 1j * gam[:nz]
print(f"  zeta'(rho) for {nz} zeros [{time.time()-t0:.0f} s]; |zeta'| at the first three: {np.abs(dz[:3]).round(3)}")
X = 10**6
mu = np.ones(X + 1, dtype=np.int8)
mu[0] = 0
for p in primes_upto(X):
    mu[p::p] *= -1
    mu[p * p::p * p] = 0
M = np.cumsum(mu)
xs = np.geomspace(100, X, 300)
real = np.array([M[int(x)] for x in xs])
pred = np.array([2 * ((x ** rho) / (rho * dz)).real.sum() - 2 for x in xs])
print(f"  corr(M(x), zeros sum) over 300 log-spaced x in [100, 1e6]: {np.corrcoef(real, pred)[0,1]:+.3f}; regression slope {np.polyfit(pred, real, 1)[0]:.3f}; rms residual / sqrt x {np.sqrt((((real - pred) / np.sqrt(xs))**2).mean()):.3f} (rms M/sqrt x {np.sqrt(((real/np.sqrt(xs))**2).mean()):.3f})")
for x in (1000, 10**4, 10**5, 10**6):
    p_ = 2 * ((x ** rho) / (rho * dz)).real.sum() - 2
    print(f"    x={x:8d}: M = {M[x]:5d}, zeros sum {p_:+8.1f}")

print("E-A7 Moebius spectrum sum mu(n) n^(-1/2-it) near the first zeros (X = 1e6) against 1/|zeta(1/2+it)|:")
n = np.arange(1, X + 1, dtype=np.float64)
w = mu[1:].astype(np.float64) / np.sqrt(n)
ln = np.log(n)
for g in gam[:4]:
    tg = np.array([g - 0.3, g - 0.1, g - 0.03, g, g + 0.03, g + 0.1, g + 0.3])
    S = np.array([(w * np.exp(-1j * t * ln)).sum() for t in tg])
    inv = np.array([1 / abs(complex(mp.zeta(mp.mpc(0.5, t)))) if t != g else float('inf') for t in tg])
    print(f"  gamma={g:.4f}: |S| at t-g = -0.3,-0.1,-0.03,0,+0.03,+0.1,+0.3: " + " ".join(f"{abs(v):6.2f}" for v in S) + "   1/|zeta|: " + " ".join(f"{v:6.2f}" if np.isfinite(v) else "  inf " for v in inv))
print(f"[{time.time()-t0:.0f} s]")
