"""Independent reference for prime17.html: Psi(x, y) by building every y-smooth number up to x (depth-first products of
primes <= y, not a sieve), x rho(u) and de Bruijn's Lambda (rho by the trapezoid rule, step 1e-4; Lambda by
Gauss-Legendre with 6 nodes on every unit interval).   usage: ref17.py x y1,y2,... out.json"""
import sys, json, math
import numpy as np
from sympy import primerange
x = int(sys.argv[1]); ys = [int(v) for v in sys.argv[2].split(',')]
h = 1e-4; U = np.arange(0, 20 + h / 2, h); rho = np.ones_like(U); k = int(round(1 / h))
for i in range(k + 1, len(U)): rho[i] = rho[i - 1] - h * (rho[i - k] / U[i] + rho[i - 1 - k] / U[i - 1]) / 2
R = lambda u: np.where(np.asarray(u) < 0, 0.0, np.interp(u, U, rho))
def psi(y):
    P = list(primerange(2, y + 1)); cnt = 1; stack = [(1, 0)]
    while stack:
        v, i = stack.pop()
        for j in range(i, len(P)):
            w = v * P[j]
            if w > x: break
            cnt += 1; stack.append((w, j))
    return cnt
GL, GW = np.polynomial.legendre.leggauss(6)
def lam(y):
    ly = math.log(y); u = math.log(x) / ly; m = np.arange(1, x + 1, dtype=float)
    s1 = float(np.sum(R(u - np.log(m) / ly) / m)); a = m[:-1]
    s2 = float(np.sum(a * sum(w * R(u - np.log(a + (g + 1) / 2) / ly) / (a + (g + 1) / 2) ** 2 for g, w in zip(GL, GW)) / 2))
    return x * (s1 - s2)
out = {}
for y in ys:
    u = math.log(x) / math.log(y); out[str(y)] = dict(ps=psi(y), r=round(x * float(R(u)), 3), L=round(lam(y), 3))
json.dump(out, open(sys.argv[3], 'w')); print(out)
