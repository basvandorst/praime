import sys, time
import numpy as np
import mpmath as mp
from math import log, sqrt, pi
sys.path.insert(0, "/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins")
from explib import zeros, primes_upto, lam_array

t0 = time.time()
gam = zeros()
X = 10**6
lam = np.ones(X + 1, dtype=np.int8); lam[0] = 0
for p in primes_upto(X):
    q = p
    while q <= X:
        lam[q::q] *= -1; q *= p
Lx = np.cumsum(lam)
print("E-A9 Liouville sums from the zeros (Perron on zeta(2s)/zeta(s)): L(x) = sqrt(x)/zeta(1/2) + sum_rho x^rho zeta(2 rho)/(rho zeta'(rho)) + ...  (1000 zeros)")
mp.mp.dps = 15
nz = 1000
rho = 0.5 + 1j * gam[:nz]
dz = np.array([complex(mp.zeta(mp.mpc(0.5, g), derivative=1)) for g in gam[:nz]])
z2 = np.array([complex(mp.zeta(mp.mpc(1.0, 2 * g))) for g in gam[:nz]])
z12 = float(mp.zeta(0.5))
print(f"  zeta(1/2) = {z12:.5f}; main term sqrt(x)/zeta(1/2) = {1/z12:.4f} sqrt x explains Polya's negativity  [{time.time()-t0:.0f} s]")
xs = np.geomspace(1000, X, 300)
real = np.array([Lx[int(x)] for x in xs])
main = np.sqrt(xs) / z12
osc = np.array([2 * ((x ** rho) * z2 / (rho * dz)).real.sum() for x in xs])
pred = main + osc
print(f"  corr(L(x), prediction) {np.corrcoef(real, pred)[0,1]:+.4f}; corr(L - main, oscillatory zeros part) {np.corrcoef(real - main, osc)[0,1]:+.4f}; regression slope of (L - main) on osc {np.polyfit(osc, real - main, 1)[0]:.3f}")
print(f"  rms (L - pred)/sqrt x = {sqrt((((real-pred)/np.sqrt(xs))**2).mean()):.4f} against rms (L - main)/sqrt x = {sqrt((((real-main)/np.sqrt(xs))**2).mean()):.4f} and rms L/sqrt x = {sqrt(((real/np.sqrt(xs))**2).mean()):.4f}")
for x in (10**4, 10**5, 10**6):
    i = np.searchsorted(xs, x) - 1
    print(f"    x={int(xs[i]):7d}: L = {real[i]:6d}, main {main[i]:8.1f}, main + zeros {pred[i]:8.1f}")

print("E-D5 Fujii's Goldbach formula at larger X (FFT to 2^24) and the scaling of the remainder:")
M = 2**24
Lam = lam_array(M)
F = np.fft.rfft(Lam, 2 * M + 2)
r = np.fft.irfft(F * F, 2 * M + 2)[:M + 1]
R = np.cumsum(r)
rho_all = 0.5 + 1j * gam
xs = np.geomspace(10**4, M, 500).astype(int)
D = np.array([R[x] - x * x / 2 for x in xs])
Z = np.array([-2 * (2 * (x ** (rho_all + 1) / (rho_all * (rho_all + 1))).real.sum()) for x in xs])
print(f"  corr {np.corrcoef(D, Z)[0,1]:+.4f}, slope {np.polyfit(Z, D, 1)[0]:.4f} over 500 X in [1e4, 1.68e7]")
res = D - Z
for a, b in ((10**4, 10**5), (10**5, 10**6), (10**6, 10**7), (10**7, M)):
    m = (xs >= a) & (xs < b)
    print(f"    X in [{a:.0e},{b:.0e}): mean (D - Z)/X = {(res[m]/xs[m]).mean():+.4f}, rms (D - Z)/X = {sqrt(((res[m]/xs[m])**2).mean()):.4f}, rms (D - Z)/X^(4/3) = {sqrt(((res[m]/xs[m]**(4/3))**2).mean()):.4f}, rms Z/X^(3/2) = {sqrt(((Z[m]/xs[m]**1.5)**2).mean()):.4f}")
print(f"[{time.time()-t0:.0f} s]")
