"""Case B: (1) does the referee's refined max miss maxima? brute force b-grid at random a;
(2) re-integrate the referee's g with the certificate's own float scheme (linear trapezoid grid) to explain
the 1.4e-6 difference between the two float J values."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from fractions import Fraction as Fr
from scipy.interpolate import CubicSpline
from ref_core import Mixture

HERE = os.path.dirname(os.path.abspath(__file__))
meas = json.load(open(os.path.join(HERE, '..', 'meas_K48.json')))
mr = [Fr(x) for x in meas['radii']]; mw = [Fr(x) for x in meas['mu']]; nwv = [Fr(x) for x in meas['nu']]
mix = Mixture([float(r) for r, w in zip(mr, mw) if w], [float(w) for w in mw if w],
              [float(r) for r, w in zip(mr, nwv) if w], [float(w) for w in nwv if w])
Psi = mix.Psi
dat = np.load(os.path.join(HERE, 'ref_g_K48.npz'))
A, G = dat['A'], dat['G']
order = np.argsort(A); A = A[order]; G = G[order]
A, idx = np.unique(A, return_index=True); G = G[idx]

rng = np.random.default_rng(99)
# (1) brute force at random a (log-uniform), compare with spline of referee g (and direct)
xs = np.exp(rng.uniform(np.log(2**-12), np.log(600), 120))
spl = CubicSpline(np.log(A), np.log(G))
worst = 0.0
for x in xs:
    bmax = x if x <= 1 else min(x, x / (x - 1))
    b = np.linspace(0, bmax, 40001)
    h = Psi(np.full_like(b, x), b) - 0.5 * np.maximum(Psi(b, b), 0)
    gb = h.max()
    gs = np.exp(spl(np.log(x)))
    rel = (gb - gs) / gb
    worst = max(worst, rel)
print(f"(1) brute force (40001 b-points) minus spline of referee g, max relative = {worst:.3e} over 120 random a")

# (2) the certificate's float scheme: trapezoid on linspace(1e-4,4,4000) + log grid to 600, with referee g
xa = np.concatenate([np.linspace(1e-4, 4, 4000), np.exp(np.linspace(np.log(4), np.log(600), 800))[1:]])
gx = np.exp(spl(np.log(xa)))
kap = 0.031770351651336
Jtheirs_scheme = 4 * np.pi * (np.trapezoid(gx / xa**5, xa) + gx[0] / xa[0]**5 * xa[0] / 2 + kap / 600)
# same g, accurate quadrature: Simpson on fine log grid
X = np.linspace(np.log(1e-4), np.log(600), 400001)
F = np.exp(spl(X)) * np.exp(-4 * X)
h = X[1] - X[0]
S = h / 3 * (F[0] + F[-1] + 4 * F[1:-1:2].sum() + 2 * F[2:-1:2].sum())
Jacc = 4 * np.pi * (S + G[0] / A[0]**6 * 1e-8 / 2 + kap / 600)
D = 263470260835153127051917226381 / 397693215151001725278407884800
c = lambda J: 1.5 * (2 * J * D * D) ** (1 / 3)
print(f"(2) referee g integrated with the certificate's float scheme: J = {Jtheirs_scheme:.8f} (c = {c(Jtheirs_scheme):.8f})")
print(f"    referee g, Simpson on 400001 log points:                  J = {Jacc:.8f} (c = {c(Jacc):.8f})")
