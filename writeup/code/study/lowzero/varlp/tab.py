"""tab.py: extend the tables of bklib.Grid / tplib (which stop at eps = 7.6) to eps in [0, EMAX] (raw units).
Same definitions as bklib.Grid and tplib/gtab.py:
  Z2 = eps^2 |zeta(1+i eps)|^2, Z2m1 = (Z2-1)/eps^2, regd = Re (zeta'/zeta)'(1+i eps) + 1/eps^2,
  g = zeta'/zeta(1+i eps) + 1/(i eps),
  Phi = prod_p (1 - w_p), w_p = (p^{i eps}-1)^2/(p-1)^2 (p <= PMAX), Phim1 = (Phi-1)/eps^2,
  d2lnPhiDiag = 2 sum_p sum_{m>=2} (m-1) ln^2 p p^-m cos(m eps ln p).
usage: tab.py EMAX H PMAX OUT.npz
"""
import sys, time, math
import numpy as np, mpmath as mp
sys.path.insert(0, 'lib')
import bklib
mp.mp.dps = 18
EMAX = float(sys.argv[1]); H = float(sys.argv[2]); PMAX = int(sys.argv[3]); OUT = sys.argv[4]
eps = np.arange(0, EMAX + H / 2, H); ev = eps.copy(); ev[0] = 1e-7
n = len(eps); Z2 = np.empty(n); Z2m1 = np.empty(n); regd = np.empty(n); g = np.empty(n, complex)
t = time.time()
for i, e in enumerate(ev):
    s = mp.mpc(1, e); z = mp.zeta(s); z1 = mp.zeta(s, derivative=1); z2 = mp.zeta(s, derivative=2); ee = mp.mpf(e)
    a = abs(1j * ee * z) ** 2
    Z2[i] = float(a); Z2m1[i] = float((a - 1) / ee ** 2)
    regd[i] = float(mp.re(z2 / z - (z1 / z) ** 2) + 1 / ee ** 2)
    g[i] = complex(z1 / z + 1 / (1j * ee))
g[0] = float(mp.euler)
print('zeta part', time.time() - t, flush=True)
P = bklib.primes_upto(PMAX).astype(float); lnP = np.log(P)
lnPhi = np.empty(n, complex); d2 = np.empty(n)
mm = []
for m in range(2, 64):
    sel = P ** m < 1e18
    if not sel.any(): break
    mm.append((m, lnP[sel], (m - 1) * lnP[sel] ** 2 * P[sel] ** (-float(m))))
CH = 100
for j in range(0, n, CH):
    e = ev[j:j + CH][:, None]
    w = -4 * np.sin(e * lnP / 2) ** 2 * np.exp(1j * e * lnP) / (P - 1) ** 2
    lnPhi[j:j + CH] = bklib.clog1p(-w).sum(axis=1)
    acc = np.zeros(e.shape[0])
    for m, lp, c in mm: acc += (c * np.cos(m * e * lp)).sum(axis=1)
    d2[j:j + CH] = 2 * acc
Phi = np.exp(lnPhi); Phim1 = bklib.cexpm1(lnPhi) / ev ** 2
np.savez(OUT, eps=eps, Z2=Z2, Z2m1=Z2m1, regd=regd, g=g, Phi=Phi, Phim1=Phim1, d2lnPhiDiag=d2, pmax=PMAX)
print('done', time.time() - t, n)
# consistency with the repo tables on [0, 7.6]
G0 = np.load('lib/grid.npz'); gt = np.load('lib/gtab.npz')
k = np.searchsorted(eps, 7.6 - 1e-9); r = G0["eps"] / H; idx = np.round(r).astype(int); ok = (idx <= n - 1) & (np.abs(r - idx) < 1e-6)
for key in ('Z2', 'Z2m1', 'regd', 'Phi', 'd2lnPhiDiag'):
    a = G0[key][ok]; b = {'Z2': Z2, 'Z2m1': Z2m1, 'regd': regd, 'Phi': Phi, 'd2lnPhiDiag': d2}[key][idx[ok]]
    print('max |diff| vs repo grid.npz', key, float(np.abs(a - b).max()))
r = gt["eps"] / H; idx = np.round(r).astype(int); ok = (idx <= n - 1) & (np.abs(r - idx) < 1e-6)
print('max |diff| vs repo gtab.npz g', float(np.abs(gt['g'][ok] - g[idx[ok]]).max()))
