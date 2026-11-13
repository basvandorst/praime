"""Primes from the logarithms of L-values: every prime below X, without zeros, without a sieve, without division.

For s > 1 and a Dirichlet character chi, log L(s, chi) = sum_m (Lambda(m)/ln m) chi(m) m^-s (principal branch: the
series has |value| <= log zeta(s) < pi). For a prime q and q not dividing n, orthogonality gives
  C(n) := (1/(q-1)) sum_{chi mod q} chibar(n) log L(s, chi) = sum_{m = n (mod q)} (Lambda(m)/ln m) m^-s.
For n < q the term m = n gives n^-s if n is prime, n^-s/k if n = p^k (k >= 2), 0 otherwise; the rest of the class
(m = n + jq) adds tail(n) in [0, n^s ((n+q)^-s + (n+q)^(1-s)/(q(s-1)))] after multiplying by n^s. So
  D(n) = n^s C(n): prime -> D >= 1, not prime -> D <= 1/2 + tailmax.
Computation: L(s, chi_k) = q^-s sum_j omega^(jk) zeta(s, g^j/q) (one FFT in discrete-log order, g a primitive root),
then log, then one more FFT gives C for every class. Nothing is crossed off and nothing is divided: the primes
come out of the logarithm (the Euler product). Usage: charlog.py X s qfactor [longdouble]
"""
import sys, math, time
import numpy as np
import scipy.fft as sfft
from scipy.special import zeta
from sympy import nextprime, primitive_root, factorint

X = int(sys.argv[1]); s = float(sys.argv[2]); qf = float(sys.argv[3]); LD = len(sys.argv) > 4
q = int(nextprime(int(qf * X)))
t0 = time.time()
g = int(primitive_root(q))
# discrete-log order pw[j] = g^j mod q, vectorised by blocks: pw[j + B] = pw[j] * g^B mod q
B = 1 << 16
pw = np.empty(q - 1, dtype=np.int64 if q > 2**31 else np.int32); v = 1
for j in range(min(B, q - 1)):
    pw[j] = v; v = v * g % q
gB = pow(g, B, q)
for start in range(B, q - 1, B):
    end = min(start + B, q - 1)
    pw[start:end] = (pw[start - B:end - B].astype(np.int64) * gB) % q
ind = np.full(q, -1, dtype=pw.dtype); ind[pw] = np.arange(q - 1, dtype=pw.dtype)
assert (ind[1:] >= 0).all()   # every nonzero residue reached once: g is a primitive root
t1 = time.time()
H = zeta(s, pw / q)                                   # Hurwitz zeta(s, a/q), a = g^j mod q
del pw
if LD:
    H = H.astype(np.longdouble)
L = sfft.ifft(H, workers=1); del H                     # L(s, chi_k) = q^-s (q-1) ifft(H)_k
L *= (q - 1) * float(q) ** (-s)
np.log(L, out=L)                                       # log L(s, chi_k), principal branch
C = sfft.fft(L, overwrite_x=True, workers=1); del L    # C[j] = class sum for n = g^j mod q
C /= (q - 1)
t2 = time.time()
n = np.arange(2, X + 1)
Cn = C[ind[n]]
D = n.astype(float) ** s * Cn.real.astype(float)
imag = float(np.max(np.abs(Cn.imag.astype(float)) * n.astype(float) ** s))
tailmax = n ** s * ((n + q) ** (-s) + (n + q) ** (1 - s) / (q * (s - 1)))
prime = D > 0.75
t3 = time.time()
print(f'X = {X}, s = {s}, q = {q} (prime, q/X = {q / X:.2f}), primitive root {g}; {"long double" if LD else "double"}')
print(f'times: discrete log {t1 - t0:.1f} s, L-values + log + FFTs {t2 - t1:.1f} s, decisions {t3 - t2:.1f} s')
print(f'largest tail bound {tailmax.max():.4f}; largest |Im D| {imag:.2e}; min D over decided primes {D[prime].min():.6f}; max D over the rest {D[~prime].max():.6f}')
# validation only (not part of the method): each n tested on its own with gmpy2 (BPSW, proved below 2^64)
import gmpy2
truth = np.fromiter((gmpy2.is_prime(int(k)) for k in n), dtype=bool, count=len(n))
wrong = np.nonzero(prime != truth)[0] + 2
print(f'primes found: {int(prime.sum())}; true count {int(truth.sum())}; wrong: {len(wrong)} {wrong[:10].tolist()}')
# exact identity: D(n) - Lambda(n)/ln(n) must lie in [0, tailmax] up to rounding
isp = np.zeros(X + 1, dtype=bool); isp[2:] = truth
lam_over_ln = truth.astype(float)
pp = np.zeros(X + 1)
for p in np.nonzero(isp)[0]:
    k, v = 2, p * p
    while v <= X:
        pp[v] = 1.0 / k; k += 1; v *= p
lam_over_ln += pp[2:]
resid = D - lam_over_ln
print(f'D(n) - Lambda(n)/ln n: min {resid.min():.3e} (must be >= 0 up to rounding), max over n of (resid - tailmax) {np.max(resid - tailmax):.3e} (must be <= 0)')
np.save(f'charlog_{X}.npy', np.stack([n, D, tailmax]))
