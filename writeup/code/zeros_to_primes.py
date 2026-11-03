"""Odlyzko's zeros: curl -sL -o zeros1 http://www.dtc.umn.edu/~odlyzko/zeta_tables/zeros1 (first 100000 zeros, 1.8 MB), kept outside the repo.
Usage: python3 -I zeros_to_primes.py <path to zeros1> [N=10000]
"""
import sys, time
import numpy as np
from math import log, pi

ZEROS = sys.argv[1]
N = int(sys.argv[2]) if len(sys.argv) > 2 else 10000
gam = np.loadtxt(ZEROS, max_rows=100000)
rho = 0.5 + 1j * gam
t0 = time.time()
xs = np.arange(1.5, N + 1.0, 1.0)
psi = np.empty(len(xs))
for i in range(0, len(xs), 400):
    x = xs[i:i + 400]
    lx = np.log(x)
    term = (np.sqrt(x)[:, None] * np.exp(1j * np.outer(lx, gam)) / rho[None, :]).sum(axis=1)
    psi[i:i + 400] = x - 2 * term.real - log(2 * pi) - 0.5 * np.log1p(-1 / x**2)
print(f"psi(x) from {len(gam)} zeros at {len(xs)} half-integers up to {N}  [{time.time()-t0:.0f} s]")
inc = np.diff(psi)
n = np.arange(2, N + 1)
s = np.ones(N + 1, dtype=bool); s[:2] = False
for i in range(2, int(N**0.5) + 1):
    if s[i]:
        s[i * i::i] = False
Lam = np.zeros(N + 1)
for p in np.flatnonzero(s):
    q = p
    while q <= N:
        Lam[q] = log(p); q *= p
Lam = Lam[2:]
err = inc - Lam
print(f"recovered Lambda(n) for 2 <= n <= {N}: max |error| = {np.abs(err).max():.3f}, rms = {np.sqrt((err**2).mean()):.3f}")
pred = inc > 0.5 * np.log(n)
truth = Lam > 0
primes_true = s[2:]
print(f"prime powers: {truth.sum()} true, {pred.sum()} predicted, missed {int((truth & ~pred).sum())}, false {int((~truth & pred).sum())}")
print(f"primes <= {N}: {primes_true.sum()}; recovered as prime powers: {int((primes_true & pred).sum())}")
worst = np.argsort(-np.abs(err))[:5]
print("largest errors at n =", ", ".join(f"{n[i]} ({err[i]:+.2f})" for i in worst))
for k in (10, 100, 1000, 10000):
    if k <= N:
        m = n <= k
        print(f"  n <= {k}: max |error| {np.abs(err[m]).max():.3f}")
