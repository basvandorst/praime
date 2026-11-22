"""Independent recomputation of G at a few maxima: plain loop n <= y, Lambda from sympy.factorint, chi built from scratch."""
import json, cmath, math
from sympy import factorint, primitive_root
K = json.load(open('consts.json')); S = json.load(open('scan_100000000.json'))
def chi_fn(q, k, m):        # chi(g^j) = e^(2 pi i k j / m), q prime
    g = primitive_root(q); dlog = {pow(g, j, q): j for j in range(m)}
    return lambda n: 0 if n % q == 0 else cmath.exp(2j * math.pi * k * dlog[n % q] / m)
cases = [('chi4', lambda n: [0, 1, 0, -1][n % 4]), ('q13 ord12 k5', chi_fn(13, 5, 12)), ('q11 ord10 k3', chi_fn(11, 3, 10)),
         ('q13 ord3 k4', chi_fn(13, 4, 12)), ('q7 ord6 k1', chi_fn(7, 1, 6)), ('q13 ord2 k6', chi_fn(13, 6, 12))]
for name, chi in cases:
    y = S[name]['at']; Y = int(y)
    if Y > 2e6: print(name, 'skip y =', y); continue
    k = K[name]; dl1 = complex(*k['dl1']); c0 = complex(*k['c0'])
    G = y * dl1
    for n in range(2, Y + 1):
        f = factorint(n)
        if len(f) == 1:
            p = next(iter(f)); G += math.log(p) * chi(n) * (y / n - 1)
    c = c0 + (0 if k['a'] else math.log(y) + 1)
    r = abs(G - c) / math.sqrt(y)
    print('%-14s y=%d ratio=%.6f scan=%.6f diff=%.1e' % (name, Y, r, S[name]['max'], abs(r - S[name]['max'])))
