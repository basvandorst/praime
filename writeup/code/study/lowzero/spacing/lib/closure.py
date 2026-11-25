"""The wave (deformation) closure: the zeros of the L-function with p a good prime are the zeros with p missing,
re-mapped by p's explicit-formula wave.  In unfolded units (mean spacing 1, L = ln(q t/2pi)):
   y = N(x) = x + dN_p(x),   dN_p(x) = -(1/pi) sum_m p^(-m/2) sin(m psi)/m = (1/pi) arg(1 - p^(-1/2) e^(i psi)),
   psi = k x + phi,  k = 2 pi ln p / L,  phi uniform (chi(p) and the height are absorbed in phi).
x: positions of the with-p zeros, y: positions of the base (p missing).  Exact in the wave (all harmonics); the closure
assumption is that the base is independent of phi.  Averages over phi by the trapezoid rule (spectrally exact)."""
import math
import numpy as np


class Wave:
    def __init__(self, p, L, M=60):
        self.p, self.L = p, L
        self.k = 2 * math.pi * math.log(p) / L
        self.a = p ** -0.5
        self.m = np.arange(1, M + 1)
        self.c = self.a ** self.m

    def dN(self, x, phi):
        psi = self.k * np.asarray(x) + phi
        return np.angle(1 - self.a * np.exp(1j * psi)) / math.pi

    def rho(self, x, phi):
        """N'(x) = 1 + d/dx dN"""
        psi = self.k * np.asarray(x) + phi
        z = self.a * np.exp(1j * psi)
        # d/dpsi arg(1 - z) = Im( -i z/(1 - z) ) = -Re(z/(1-z))
        return 1 - (self.k / math.pi) * np.real(z / (1 - z))

    def N(self, x, phi):
        return np.asarray(x) + self.dN(x, phi)

    def inv_right(self, x0, u, phi, it=60):
        """s >= 0 with N(x0 + s) - N(x0) = u (Newton; N is increasing since |dN'| < 1 here)"""
        y0 = self.N(x0, phi); s = np.array(u, float, copy=True)
        for _ in range(it):
            f = self.N(x0 + s, phi) - y0 - u
            s = s - f / self.rho(x0 + s, phi)
        return s

    def inv_left(self, x0, u, phi, it=60):
        """s >= 0 with N(x0) - N(x0 - s) = u"""
        y0 = self.N(x0, phi); s = np.array(u, float, copy=True)
        for _ in range(it):
            f = y0 - self.N(x0 - s, phi) - u
            s = s - f / self.rho(x0 - s, phi)
        return s


def phases(n=64):
    return 2 * math.pi * np.arange(n) / n


def delta_count_var(base_sig2, wave, ells, nphi=128):
    """Number variance of the re-mapped process: E[Sigma2_base(l + D)] + Var(D), D = dN(x0+l) - dN(x0), x0 uniform.
    base_sig2: callable l -> Sigma2_base(l). Returns Sigma2_with(l) - Sigma2_base(l)."""
    ph = phases(nphi); out = []
    for l in ells:
        D = wave.dN(l, ph) - wave.dN(0.0, ph)
        out.append(np.mean(base_sig2(l + D)) + np.mean(D ** 2) - base_sig2(np.array([l]))[0])
    return np.array(out)
