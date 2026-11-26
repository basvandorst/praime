"""Validate wform.py on a smooth even test function f(u) = (1-(2u/L)^2)^4 on [-L/2, L/2], L = 3:
(1) Legendre-matrix value c^T D_primes c,
(2) Fourier-side value pole + (1/pi) int_0^inf |f^|^2 a(t) dt - comb  (independent of the matrix code),
(3) Weil explicit formula: 2 sum_{gamma>0, first 100 zeros} |f^(gamma)|^2.
Observed (Nov 2026): (1) 1.48721409e-08, (3) 1.48721044e-08, i.e. agreement to 4e-14 absolute
after a cancellation from size ~7."""
import numpy as np, mpmath as mp
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
import wform as W

L, N = 3.0, 40
f = lambda u: (1 - (2 * u / L) ** 2) ** 4
t, w = leggauss(3000); x = L / 2 * t; wx = w * L / 2
c = (W.phi(N, L, x) * (wx * f(x))[:, None]).sum(0)
P = W.primes_upto(np.exp(L) + 1)
v1 = c @ W.form_D(N, L, P) @ c

fhat = lambda tt: np.sum(wx * f(x) * np.cos(tt * x))
a = lambda tt: float(mp.re(mp.digamma(0.25 + 0.5j * tt))) - np.log(np.pi)
pole = 2 * np.sum(wx * f(x) * np.exp(x / 2)) * np.sum(wx * f(x) * np.exp(-x / 2))
arch = quad(lambda tt: fhat(tt) ** 2 * a(tt), 0, 300, limit=3000)[0] / np.pi
def gcorr(u):
    lo, hi = u - L / 2, L / 2
    tt, ww = leggauss(200)
    xx = (hi - lo) / 2 * tt + (hi + lo) / 2
    return np.sum(ww * (hi - lo) / 2 * f(xx) * f(xx - u))
comb = sum(wt * 2 * gcorr(pos) for pos, wt in W.comb_atoms(P, L))
v2 = pole + arch - comb
mp.mp.dps = 15
zs = [float(mp.im(mp.zetazero(n))) for n in range(1, 101)]
v3 = 2 * sum(fhat(z) ** 2 for z in zs)
print('(1) Legendre matrix :', v1)
print('(2) Fourier side    :', v2, '  [pole', pole, 'arch', arch, 'comb', comb, ']')
print('(3) zeros (100)     :', v3)
