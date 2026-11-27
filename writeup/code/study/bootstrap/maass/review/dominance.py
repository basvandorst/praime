"""Check of the parity-dominance identity: Q_odd - Q_even = (1/2) int |phihat|^2 sum_rho sech(pi (t - rho)) dt >= 0
(rho over the imaginary parts of the shifts), on the matrices: lambda_min(A_k^{eps=1} - A_k^{eps=0}) >= 0 and the
t-space value for the smooth phi of check_norm.py."""
import sys; sys.path.insert(0, "/home/bas/praime/writeup/code/study/bootstrap/maass/review")
import numpy as np, mpmath as mp, weilmat as W
mp.mp.dps = 20
geo = W.Geo(8.0, 32)
for r in (9.5337, 13.78, 15.0):
    for k in (1, 3):
        D = W.arch(geo, W.gammas(k, r, 1)) - W.arch(geo, W.gammas(k, r, 0))
        ev = np.linalg.eigvalsh(D)
        print("r=%.4f k=%d  eig(D) min %.3e max %.3e" % (r, k, ev[0], ev[-1]))
# t-space identity on one function: Legendre P_0 normalised: phihat = sqrt(1/(2a)) 2 sin(a t)/t
a = geo.a
ph = lambda t: (mp.sqrt(1 / (2 * a)) * 2 * mp.sin(a * t) / t) if t != 0 else mp.sqrt(1 / (2 * a)) * 2 * a
r = 9.5337
lhs = (W.arch(geo, W.gammas(1, r, 1)) - W.arch(geo, W.gammas(1, r, 0)))[0, 0]
rhs = 0.5 * mp.quad(lambda t: ph(t) ** 2 * (mp.sech(mp.pi * (t - r)) + mp.sech(mp.pi * (t + r))), [-mp.inf, -r, 0, r, mp.inf])
print("P_0 entry: matrix %.12e  sech formula %.12e" % (lhs, float(rhs)))
