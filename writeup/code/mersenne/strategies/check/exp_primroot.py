"""Exploratory: least primitive root of the Mersenne exponents vs all primes (Artin). One test, labelled as such."""
import sympy, math
E = [3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243,110503,132049,216091,756839,859433,1257787,1398269,2976221,3021377,6972593,13466917,20996011,24036583,25964951,30402457,32582657,37156667,42643801,43112609,57885161,74207281,77232917,82589933,136279841]
g = {p: sympy.primitive_root(p, smallest=True) for p in E}
k2 = sum(1 for p in E if g[p] == 2)
# null: 2 is a primitive root mod p needs p = +-3 mod 8 (half the primes) and then about 2A of those; condition on each exponent's residue mod 8
A = 0.3739558136
m3 = [p for p in E if p % 8 in (3, 5)]
k2m = sum(1 for p in m3 if g[p] == 2)
from scipy.stats import binomtest
print('exponents 3..136279841:', len(E), '; least primitive root 2:', k2, f'({100*k2/len(E):.1f}%, Artin {100*A:.1f}%)', '; binomial p', round(binomtest(k2, len(E), A).pvalue, 3))
print('exponents = +-3 mod 8:', len(m3), '; with 2 as primitive root:', k2m, f'(expected share about 2A = {200*A:.1f}%)', '; binomial p', round(binomtest(k2m, len(m3), 2 * A).pvalue, 3))
print('least primitive roots:', [g[p] for p in E])
