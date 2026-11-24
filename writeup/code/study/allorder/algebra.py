"""Exact-algebra checks (no primes needed beyond the pooled matrices of the stretches):
 (1) double centring of the within marks z_d(a, b) = B_d,ab / sqrt(pi_a P_ab): sum_b P_ab z_d(a,b) = 0 and sum_a pi_a P_ab z_d(a,b) = 0;
 (2) hence every class-level four-pair term vanishes: for a pairing that couples the class of one prime of event 1 with the
     class of one prime of event 2 through a weight w(c, c') (any function), E[z1 z2 w] = 0 when one of the two marks is
     averaged over its free class; checked for random w;
 (3) parity of the character-pair directions under (a, b) -> (-a, -b) and the number of odd directions for q = 3..13;
 (4) Jacobi sums J(chi_k, chi_l) = sum_t chi_k(t) conj chi_l(1 + t) mod 5 (magnitudes 1 or sqrt 5)."""
import numpy as np, json, math, cmath, sys
sys.path.insert(0, "./")
from common import *
L = json.load(open('lagcov.json'))
rng = np.random.default_rng(0)
for key in ('B_5', 'B_7', 'A_5'):
    E = L[key]; q = int(key.split('_')[1]); P = np.array(E['P']); pi = np.array(E['pi']); C = len(pi)
    B, lab, par = basis(q, P, pi)
    Z = B.reshape(C, C, -1) / np.sqrt(pi[:, None, None] * P[:, :, None])
    row = np.abs(np.einsum('ab,abd->ad', P, Z)).max(); col = np.abs(np.einsum('a,ab,abd->bd', pi, P, Z)).max()
    # pairing check: E over (a,b) x (c,d) of z1(a,b) z2(c,d) w(x, y) with x in {a,b}, y in {c,d}
    J = pi[:, None] * P
    worst = 0
    for _ in range(5):
        w = rng.normal(size=(C, C))
        for i1 in (0, 1):
            for i2 in (0, 1):
                for d in range(Z.shape[2]):
                    s = 0.0
                    for a in range(C):
                        for b in range(C):
                            x = (a, b)[i1]
                            for c in range(C):
                                for dd in range(C):
                                    y = (c, dd)[i2]
                                    s += J[a, b] * J[c, dd] * Z[a, b, d] * Z[c, dd, d] * w[x, y]
                    worst = max(worst, abs(s))
    print('%s: max |E[z|first class]| %.1e  max |E[z|second class]| %.1e  max |class-level pairing term| %.1e' % (key, row, col, worst))
for q in (3, 5, 7, 11, 13):
    ch = chars(q); C = len(red(q)); im1 = red(q).index(q - 1)
    sgn = [int(round(ch[k][im1].real)) for k in range(C)]
    odd = sum(1 for k in range(1, C) for l in range(1, C) if sgn[k] * sgn[l] == -1)
    print('q=%2d: chi_k(-1) = %s ; odd directions (chi_k(-1) chi_l(-1) = -1): %d of %d' % (q, sgn[1:], odd, (C - 1) ** 2))
ch = chars(5); r = red(5)
for k in range(1, 4):
    print('  q=5 Jacobi |J(chi_%d, chi_l)|, l=1..3:' % k, ['%.3f' % abs(sum(ch[k][r.index(t)] * np.conj(ch[l][r.index((1 + t) % 5)]) for t in r if (1 + t) % 5)) for l in range(1, 4)])
