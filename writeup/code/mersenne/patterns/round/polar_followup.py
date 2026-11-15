"""Polar layout: (a) which exponents drive the angle concentration (post hoc, by size threshold);
(b) the rendered bitmap: occupancy and prime-pixel fraction against radius (explains the disc)."""
import numpy as np, common as C
isp = C.sieve(); P, w = C.null_weights(isp)
null = C.sample_null(P, w, 3000, 31, seed=99)
E = C.MEXP
def ray(v):
    return np.abs(np.exp(1j*v.astype(float)).sum())**2/len(v)
print('angles (rad mod 2pi, in (-pi,pi]):', {int(k): round(float(np.angle(np.exp(1j*k))), 2) for k in E[:14]})
print('non-exponent primes <= 31:', {k: round(float(np.angle(np.exp(1j*k))), 2) for k in (11, 23, 29)})
for T in (1, 32, 128, 1000):
    o = E[E >= T]; obs = ray(o)
    nl = np.array([ray(s[s >= T]) for s in null])
    print(f'exponents >= {T}: m={len(o)}  Rayleigh k=1 = {obs:.2f}  null mean {nl.mean():.2f}  p = {(1+np.sum(nl >= obs))/(1+len(nl)):.3f}')
# rendered bitmap: pixel keeps the largest n (page: grid[...] = n in increasing n)
pos = np.fromfile(C.DATA+'/polar.pos.bin', '<u2').reshape(-1, 2)[1:C.N+1].astype(np.int64)
print('mean angle of exponents (deg):', round(float(np.degrees(np.angle(np.exp(1j*E.astype(float)).sum()))), 1))
W = 1404; key = pos[:, 1]*W + pos[:, 0]
last = np.zeros(W*W, np.int64); last[key] = np.arange(1, C.N+1)   # numpy fancy assignment: last write wins
occ = last > 0; prm = isp[last] & occ
yy, xx = np.divmod(np.arange(W*W), W); r = np.hypot(xx+0.5-702, yy+0.5-702)
for r0 in range(0, 700, 50):
    s = (r >= r0) & (r < r0+50)
    sp = 2*np.pi*(r0+25)/710
    print(f'r {r0:3d}-{r0+50:3d}: occupied {occ[s].mean():.2f}  prime pixels {prm[s].mean():.3f}  spacing of 710 spokes {sp:.2f}px')
