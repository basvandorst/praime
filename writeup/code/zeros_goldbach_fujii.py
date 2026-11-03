import sys, time
import numpy as np
from math import log, pi

gam = np.loadtxt(sys.argv[1], max_rows=100000)
rho = 0.5 + 1j * gam
M = 2**21
t0 = time.time()
s = np.ones(M + 1, dtype=bool); s[:2] = False
for i in range(2, int(M**0.5) + 1):
    if s[i]: s[i * i::i] = False
Lam = np.zeros(M + 1)
for p in np.flatnonzero(s):
    q = p
    while q <= M:
        Lam[q] = log(p); q *= p
F = np.fft.rfft(Lam, 2 * M + 2)
r = np.fft.irfft(F * F, 2 * M + 2)[:M + 1]
R = np.cumsum(r)
print(f"r(N) = sum Lambda(m)Lambda(N-m) by FFT up to {M}  [{time.time()-t0:.0f} s]")
print("Fujii: sum_{N<=X} r(N) = X^2/2 - 2 sum_rho X^(rho+1)/(rho(rho+1)) + lower order")
print("     X        D(X) = sum r - X^2/2     Z(X) = -2 Re sum X^(rho+1)/(rho(rho+1))   D/Z    (D-Z)/X^(4/3)")
Xs = [10**4, 2 * 10**4, 5 * 10**4, 10**5, 2 * 10**5, 5 * 10**5, 10**6, 2 * 10**6]
D, Z = [], []
for X in Xs:
    d = R[X] - X * X / 2
    zz = -2 * (2 * (X ** (rho + 1) / (rho * (rho + 1))).real.sum())
    D.append(d); Z.append(zz)
    print(f"  {X:8d}   {d:+16.0f}   {zz:+16.0f}   {d/zz:6.3f}   {(d-zz)/X**(4/3):+8.3f}")
xs = np.geomspace(10**4, 2 * 10**6, 400).astype(int)
Dd = np.array([R[x] - x * x / 2 for x in xs])
Zz = np.array([-2 * (2 * (x ** (rho + 1) / (rho * (rho + 1))).real.sum()) for x in xs])
print(f"over 400 log-spaced X in [1e4, 2e6]: corr(D, Z) = {np.corrcoef(Dd, Zz)[0,1]:+.4f}, regression slope {np.polyfit(Zz, Dd, 1)[0]:.3f}, rms (D - Z)/X^(4/3) = {np.sqrt((((Dd - Zz) / xs**(4/3))**2).mean()):.3f}")
print(f"with the factor 1 instead of 2 the slope would read {np.polyfit(Zz/2, Dd, 1)[0]:.3f}")
