"""Independent check of infeasibility certificates.
For a system P and level L, take the lowest eigenvector of the Legendre matrix (wform.py),
form the explicit polynomial f on [-L/2, L/2], and evaluate D_P(f * f~) again by a different route:
  pole  = 2 (int f e^{u/2}) (int f e^{-u/2})              (f real)
  arch  = (1/pi) int_0^inf |f^(t)|^2 a(t) dt,  a(t) = Re psi(1/4+it/2) - log pi   (Fourier side)
  comb  = sum_atoms w * 2 * (f*f~)(position)               (direct correlation by quadrature)
A negative value with |error| << |value| certifies that P is not feasible at level L."""
import numpy as np, mpmath as mp, sys
from numpy.polynomial.legendre import leggauss
from scipy.special import digamma
import wform as W

def a_fun(t):
    # Re psi(1/4 + i t/2) - log pi, vectorised through mpmath for moderate t, asymptotic for large t
    out = np.empty_like(t)
    for i, tt in enumerate(t):
        out[i] = float(mp.re(mp.digamma(mp.mpf(1) / 4 + 0.5j * tt))) - np.log(np.pi)
    return out

def certify(P, L, N=40, Tmax=3000.0):
    M = W.form_D(N, L, P)
    ev, V = np.linalg.eigh((M + M.T) / 2)
    c = V[:, 0]
    tq, wq = leggauss(4000)
    x = L / 2 * tq; wx = wq * L / 2
    fx = W.phi(N, L, x) @ c
    norm2 = np.sum(wx * fx ** 2)
    pole = 2 * np.sum(wx * fx * np.exp(x / 2)) * np.sum(wx * fx * np.exp(-x / 2))
    # Fourier side: t grid
    tg, wg = leggauss(400)
    edges = np.linspace(0, Tmax, 301)
    arch = 0.0
    for a, b in zip(edges[:-1], edges[1:]):
        t = (b - a) / 2 * tg + (a + b) / 2
        wt = wg * (b - a) / 2
        F = (np.cos(np.outer(t, x)) @ (wx * fx)) ** 2 + (np.sin(np.outer(t, x)) @ (wx * fx)) ** 2
        arch += np.sum(wt * F * a_fun(t)) / np.pi
    fp, fm = (W.phi(N, L, np.array([L / 2])) @ c)[0], (W.phi(N, L, np.array([-L / 2])) @ c)[0]
    tail = (fp ** 2 + fm ** 2) * (np.log(Tmax / (2 * np.pi)) + 1) / Tmax / np.pi
    arch += tail
    comb = 0.0
    for pos, w in W.comb_atoms(P, L):
        lo, hi = pos - L / 2, L / 2
        tt, ww = leggauss(400)
        xx = (hi - lo) / 2 * tt + (hi + lo) / 2
        fa = W.phi(N, L, xx) @ c; fb = W.phi(N, L, xx - pos) @ c
        comb += w * 2 * np.sum(ww * (hi - lo) / 2 * fa * fb)
    val = pole + arch - comb
    return ev[0], val / norm2, tail / norm2

cases = [([], 0.745, 'empty system, level 0.745 (e^L=2.1066)'),
         ([2.0], np.log(3.0488) + 0.01, '{2}, level log(3.0488)+0.01'),
         ([2.0, 3.0], np.log(5.0316) + 0.01, '{2,3}, level log(5.0316)+0.01'),
         ([], np.log(2.0), 'empty system, level log 2 (expect >0)')]
for P, L, name in cases:
    lam, val, tail = certify(P, L)
    print(f'{name}: matrix lam_min={lam:.6e}  independent Rayleigh quotient={val:.6e}  (tail part {tail:.1e})')
