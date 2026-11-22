# Independent recheck (mpmath, exact incomplete gamma) of heuristic exclusions found by region2.py
import math, numpy as np, mpmath as mp
import region2 as R
mp.mp.dps = 40
env = np.load('env1e8.npy'); W = 0.05; b = mp.log(10**8); C = mp.mpf(2) + mp.euler - mp.log(4*mp.pi)
for beta, g0, mode in [(0.9, 100.0, 'full'), (0.75, 31.0, 'full'), (0.9, 86.0, 'win')]:
    lam, nu = R.excluded(beta, g0, 'full' if mode == 'full' else (math.log(1e6), float(env[int(math.log(1e6)/W):].max())), False)
    lam = mp.mpf(lam); sig = 1+lam; Rr = sig-mp.mpf(1)/2; d0 = sig-beta
    if mode == 'full':
        J = sum(mp.mpf(max(float(c), 0))*mp.gammainc(nu, Rr*j*W, Rr*min(mp.mpf(j+1)*W, b), regularized=True) for j, c in enumerate(env))
        data = 2*J*Rr**(-nu); low_tail = 0
    else:
        a = mp.log(10**6); C0 = mp.mpf(float(env[int(math.log(1e6)/W):].max()))
        data = 2*C0*mp.gammainc(nu, Rr*a, Rr*b, regularized=True)*Rr**(-nu)
        low_tail = 2*R.S1*mp.gammainc(nu, 0, lam*a, regularized=True)*lam**(-nu)
    hi_tail = 2*R.S1*mp.gammainc(nu, lam*b, mp.inf, regularized=True)*lam**(-nu)
    realpt = R.S1*(lam**2+14**2)**(-mp.mpf(nu)/2)
    U = data + low_tail + hi_tail + realpt
    rho = mp.mpc(beta, g0); r0 = 1/abs(rho*(1-rho))
    L = r0*d0**(-nu) - (C+3*r0)*Rr**(-nu)
    print(f"beta={beta} g0={g0} {mode}: sigma={float(sig)} nu={nu}: lower {mp.nstr(L,6)} > upper {mp.nstr(U,6)} ? {L > U}  [data {mp.nstr(data,3)}, tail>Y1 {mp.nstr(hi_tail,3)}, tail<Y0 {mp.nstr(low_tail,3)}, A(sigma) {mp.nstr(realpt,3)}]")
