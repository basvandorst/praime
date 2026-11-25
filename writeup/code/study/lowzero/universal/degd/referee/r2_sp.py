"""Referee check 2: one-term stationary phase I = int H(t) X^{-1}(1/2+it) n^{-it} dt.
Phase computed directly as X^{-1}(1/2+it) (loggamma), not via the note's theta.
Claim: I = 2 pi A_d (n/Q)^{1/(2d)} e(-d (n/Q)^{1/d}) [H (1 + i kappa/t) + (i/2d) H' + (i t/2d) H''](t*) + O(t^{1/2} t^{-2}),
A_d = e^{i(phi0+pi/4)}/sqrt d, phi0 = (pi/2) sum Re mu - d pi/4, kappa = (d^2-1)/(24d) + (1/2) sum Re(mu - mu^2)   (b = sum Im mu = 0).
Window 1 (plateau): H = exp(-((t - tc)/(tc/2))^8): H' = H'' = 0 at tc, so (I/lead - 1) t -> i kappa.
Window 2 (Gaussian, off-centre): H = exp(-(t - 1.15 tc)^2 / (2 (0.25 tc)^2)): tests the H', H'' terms."""
import math, cmath
import numpy as np
from scipy.special import loggamma

def lnGR(z): return -(z / 2) * math.log(math.pi) + loggamma(z / 2)
def Xinv_line(t, Q, mu):
    s = 0.5 + 1j * t
    v = 1j * t * math.log(Q) + 0j
    for m in mu:
        v = v + lnGR(s + m) - lnGR(1 - s + np.conj(m))
    return np.exp(v)

cases = [('zeta-like d=1 mu=0, Q=1', 1.0, [0.0]), ('d=1 odd char Q=7', 7.0, [1.0]),
         ('d=2 weight-2 newform Q=11', 11.0, [0.5, 1.5]), ('d=2 Delta Q=1', 1.0, [5.5, 6.5]),
         ('d=3 Sym2 Q=121', 121.0, [1.0, 1.0, 2.0]),
         ('d=3 Maass-type mu=i(2,-.5,-1.5) Q=10', 10.0, [2j, -0.5j, -1.5j]),
         ('d=4 mu=(.5,.5,1.5,1.5) Q=249', 249.0, [0.5, 0.5, 1.5, 1.5]),
         ('d=5 mu=(0,1,1,2,3) Q=3', 3.0, [0.0, 1.0, 1.0, 2.0, 3.0])]
for name, Q, mu in cases:
    d = len(mu); muc = [complex(m) for m in mu]
    b = sum(m.imag for m in muc); assert abs(b) < 1e-12
    phi0 = math.pi / 2 * sum(m.real for m in muc) - d * math.pi / 4
    kap = (d * d - 1) / (24 * d) + 0.5 * sum((m - m * m).real for m in muc)
    A = cmath.exp(1j * (phi0 + math.pi / 4)) / math.sqrt(d)
    out1, out2 = [], []
    for tc in (100.0, 200.0, 400.0, 800.0, 1600.0):
        h = 0.02
        t = np.arange(0.02 * tc, 2.3 * tc, h)
        n = Q * (tc / (2 * math.pi)) ** d
        lead0 = 2 * math.pi * A * (n / Q) ** (1 / (2 * d)) * cmath.exp(-2j * math.pi * d * (n / Q) ** (1 / d))
        E = Xinv_line(t, Q, muc) * np.exp(-1j * t * math.log(n))
        # window 1
        H1 = np.exp(-((t - tc) / (tc / 2)) ** 8)
        I1 = h * np.sum(H1 * E)
        r1 = (I1 / lead0 - 1) * tc  # -> i kappa
        # window 2 with derivatives at tc
        s2 = 0.25 * tc; c2 = 1.15 * tc
        H2 = np.exp(-(t - c2) ** 2 / (2 * s2 * s2))
        I2 = h * np.sum(H2 * E)
        g = math.exp(-(tc - c2) ** 2 / (2 * s2 * s2)); u = (tc - c2) / s2 ** 2
        H, Hp, Hpp = g, -u * g, (u * u - 1 / s2 ** 2) * g
        first = lead0 * (H * (1 + 1j * kap / tc) + 1j / (2 * d) * Hp + 1j * tc / (2 * d) * Hpp)
        leadonly = lead0 * H
        out1.append((tc, r1)); out2.append((tc, abs(I2 / leadonly - 1) * tc, abs(I2 / first - 1) * tc * tc))
    print(f'{name}: kappa(formula) = {kap:.6f}')
    print('   plateau window, (I/lead - 1)*t at t = 100..1600: ' + ', '.join(f'{z.real:+.5f}{z.imag:+.5f}i' for _, z in out1))
    print('   Gaussian window: |I/lead-1|*t = ' + ', '.join(f'{a:.4f}' for _, a, _ in out2) +
          '; |I/first-1|*t^2 = ' + ', '.join(f'{c:.3f}' for _, _, c in out2))
    # Richardson on the last two: i kappa estimate
    z1, z2 = out1[-2][1], out1[-1][1]
    est = 2 * z2 - z1
    print(f'   kappa from Richardson (1/t extrapolation): {est.imag:.6f} (real part {est.real:.2e}); formula {kap:.6f}; diff {abs(est.imag - kap):.2e}')
