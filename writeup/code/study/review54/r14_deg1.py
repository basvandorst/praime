# Degree 1: exact finite-T prime-side value S(T) = -i^a q^{-1/2} sum_{n <= qT/2pi} Lambda(n) chi(n) e(-n/q) against the measured sums.
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+"/")
import json, math, cmath, numpy as np, sympy
from r1_sums import theta, gauss_eps, D
pred = json.load(open(D + 'predictions.json'))['Lfunctions']
for key in ('X5', 'X5d', 'X7', 'X7d'):
    P = pred[key]; q = int(P['Q']); T = float(P['T2']); epsF = complex(*P['eps'])
    order = 4 if q == 5 else 3
    best = None
    for j in range(1, order):
        if math.gcd(j, order) != 1: continue
        chi, e, a = gauss_eps(q, cmath.exp(2j * math.pi * j / order))
        if best is None or abs(e - epsF) < abs(best[1] - epsF): best = (chi, e, a)
    chi, e, a = best
    x = q * T / (2 * math.pi)
    S = 0j
    for n in range(2, int(x) + 1):
        f = sympy.factorint(n)
        if len(f) != 1: continue
        p = list(f)[0]
        if n % q == 0: continue
        S += math.log(p) * chi[n % q] * cmath.exp(-2j * math.pi * n / q)
    S = -(1j ** a) / math.sqrt(q) * S
    g = np.load(D + 'zeros/zeros_%s.npz' % key)['g']; g = np.sort(g[(g > 0) & (g <= T)])
    G = np.exp(2j * theta(g, float(q), [float(a)])).sum()
    asym = -e * (q / (q - 1)) * T / (2 * math.pi)
    print('%s: measured %.2f%+.2fi | exact prime side %.2f%+.2fi | asymptotic -eps q/phi(q) T/2pi %.2f%+.2fi | measured/exact %.3f%+.3fi ; exact/asymptotic %.3f%+.3fi'
          % (key, G.real, G.imag, S.real, S.imag, asym.real, asym.imag, (G / S).real, (G / S).imag, (S / asym).real, (S / asym).imag))
