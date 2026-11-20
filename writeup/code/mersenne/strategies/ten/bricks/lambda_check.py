"""Second computation of de Bruijn's Lambda(x, y) = x int rho(u - v) d(floor(y^v)/y^v): a Stieltjes sum on a grid of
2*10^7 + 1 points in v (rho at the midpoints, jumps of floor(t)/t taken from the grid), for x = 10^6, y = 100, 1000, 10000;
and the sanity check rho = 1 -> Lambda = x.   usage: lambda_check.py"""
import math
import numpy as np
h = 1e-4; U = np.arange(0, 20 + h / 2, h); rho = np.ones_like(U); k1 = int(round(1 / h))
for i in range(k1 + 1, len(U)): rho[i] = rho[i - 1] - h * (rho[i - k1] / U[i] + rho[i - 1 - k1] / U[i - 1]) / 2
R = lambda u: np.where(np.asarray(u) < 0, 0.0, np.interp(u, U, rho))
one = lambda u: np.where(np.asarray(u) < 0, 0.0, 1.0)
def grid(x, y, Rf, K=20_000_001):
    ly = math.log(y); u = math.log(x) / ly; v = np.linspace(0, u, K); t = np.exp(v * ly); F = np.floor(t + 1e-9) / t
    return x * (float(np.sum(Rf(u - (v[1:] + v[:-1]) / 2) * np.diff(F))) + float(Rf(u)) * F[0])
x = 10 ** 6
print('rho = 1: Lambda / x = %.7f' % (grid(x, 1000, one) / x))
for y in (100, 1000, 10000): print('x = 10^6, y = %5d: Lambda = %.1f' % (y, grid(x, y, R)))
