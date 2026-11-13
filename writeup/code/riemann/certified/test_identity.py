"""Step 1: check the closed forms and the exact explicit-formula identity on real zeros."""
import sys, math
import os
import numpy as np
from scipy.integrate import quad
from kb import K, Khat, P, K0
from sympy import factorint

S = os.environ.get('WORK', os.getcwd()) + '/'
G = np.fromfile(S + 'z2m.f64')

def lam(m):
    f = factorint(m)
    return math.log(next(iter(f))) if len(f) == 1 else 0.0

a, beta = 1.9, 10.0
# closed form of Khat against quadrature
for w in [0.0, 1.0, 3.0, 5.2, 5.3, 6.0, 10.0, 40.0]:
    num = quad(lambda v: float(K(v, a, beta)) * math.cos(w * v), -a, a, limit=400)[0]
    print(f'Khat({w}) closed {float(Khat(w, a, beta)):.10e}  quad {num:.10e}')
for x in [100.0, 1009.0]:
    num = quad(lambda v: float(K(v, a, beta)) * math.exp(v / (2 * x)), -a, a, limit=400)[0]
    print(f'P({x}) closed {float(P(x, a, beta)):.12e} quad {num:.12e}')

# identity at several x (real numbers too)
for (x, a, beta) in [(1009.0, 1.9, 10.0), (1010.0, 1.9, 10.0), (1009.37, 1.9, 10.0), (10007.0, 1.9, 12.0), (5000.5, 30.0, 14.0), (99991.0, 1.9, 14.0)]:
    lo, hi = x * math.exp(-a / x), x * math.exp(a / x)
    ps = sum(lam(m) * math.sqrt(x / m) * float(K(x * math.log(m / x), a, beta)) for m in range(max(2, int(lo)), int(hi) + 2))
    w = G / x
    zs_terms = np.cos(G * math.log(x)) * Khat(w, a, beta)
    zs = float(P(x, a, beta) - (2 / math.sqrt(x)) * np.sum(zs_terms))
    T = beta / a * x
    print(f'x={x} a={a} beta={beta}: prime side {ps:.9f}  zero side {zs:.9f}  diff {zs-ps:.2e}  (K(0)={K0(beta):.4e}; band edge T = beta x/a = {T:.0f}, zeros used up to {G[-1]:.0f})')
