"""Lower bound for local linear prime tests under the random-phase law.
For a kernel f on [-1, 1] (normalized: u = v/a) with f(0) = 1, the real error is E^2 = (n c/(q-1)) (L_w/pi) int_s^inf |fhat|^2,
L_w ~ ln(q s/(2 pi c)). The smallest out-of-band energy R*(s) = min int_{|sigma|>s} |fhat|^2 over f with f(0) = 1 is
1/(e0^T M^-1 e0) for the quadratic form M = 2 pi I - B_s, (B_s)_ij = h * 2 sin(s(u_i - u_j))/(u_i - u_j) (Parseval).
Compared with the Kaiser-Bessel windows used (best beta for each s), and turned into the smallest s at which the worst
of N = 10^4 numbers stays below (1/4) ln n."""
import math
import numpy as np
from scipy.special import iv
from scipy.integrate import quad
from kbm import Z as Zkb
n = 10 ** 4; Nn = 10 ** 4; margin = 0.25 * math.log(n); m = 4
h = 2 / 600; u = np.arange(-1 + h / 2, 1, h)
D = u[:, None] - u[None, :]
def Rstar(s):
    with np.errstate(invalid='ignore', divide='ignore'):
        B = np.where(D == 0, 2 * s, 2 * np.sin(s * D) / np.where(D == 0, 1, D))
    M = h * (2 * math.pi * np.eye(len(u)) - h * B)          # int |fhat|^2 out of band = f^T M f (f sampled, f(0) ~ value at centre)
    # constraint f(0) = 1: interpolate f(0) by the two centre samples
    e0 = np.zeros(len(u)); i0 = len(u) // 2; e0[i0 - 1] = e0[i0] = 0.5
    y = np.linalg.solve(M, e0)
    return 1.0 / (e0 @ y)
def R_KB(s, beta):
    C0 = math.sqrt(2 * math.pi) * beta ** m; K0 = iv(m, beta)
    return 2 * quad(lambda x: (C0 * float(Zkb(x, beta, m)) / K0) ** 2, s, s + 500, limit=500)[0]
Lw = math.log(n * 5 / (2 * math.pi))
print('s     optimal out-of-band energy R*   best Kaiser-Bessel (beta)        worst error / margin: optimal   Kaiser-Bessel')
for s in (2.0, 2.5, 3.0, 3.6, 4.0, 4.8, 6.0, 7.2, 9.6):
    Ro = Rstar(s)
    best = min((R_KB(s, b), b) for b in np.arange(0.5, s / 1.2 + 0.01, 0.25))
    Eo = math.sqrt(Lw / (2 * math.pi) * Ro) * math.sqrt(2 * math.log(Nn)) / margin
    Ek = math.sqrt(Lw / (2 * math.pi) * best[0]) * math.sqrt(2 * math.log(Nn)) / margin
    print(f'{s:4.1f}  {Ro:.3e}                      {best[0]:.3e} ({best[1]:.2f})              {Eo:.3f}            {Ek:.3f}')
