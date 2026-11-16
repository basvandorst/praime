"""The Mersenne exponents by class mod 12, against (a) the plain heuristic (Wagstaff mod 4), (b) theory 1 (lanes from the small
factors), and the size of Chebyshev's class-1 lag among all primes in the same range."""
import numpy as np, math
E = [5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253, 4423, 9689, 9941, 11213, 19937,
     21701, 23209, 44497, 86243, 110503, 132049, 216091, 756839, 859433, 1257787, 1398269, 2976221, 3021377, 6972593,
     13466917, 20996011, 24036583, 25964951, 30402457, 32582657, 37156667, 42643801, 43112609, 57885161, 74207281,
     77232917, 82589933, 136279841]
B = 136279842
G = math.exp(0.5772156649015329)
s = np.ones(B // 2, bool); s[0] = False
for i in range(3, int(B ** 0.5) + 1, 2):
    if s[i // 2]: s[i * i // 2::i] = False
P = 2 * np.nonzero(s)[0] + 1
P = P[P > 3]
LANES = [1, 5, 7, 11]
sizes = np.log([1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8])
LF = np.array([[1.153, 1.038, 1.203, 0.606], [1.100, 1.020, 1.132, 0.747], [1.074, 1.013, 1.098, 0.815], [1.059, 1.009, 1.078, 0.854],
               [1.049, 1.007, 1.065, 0.879], [1.042, 1.006, 1.055, 0.897], [1.036, 1.005, 1.048, 0.911]])
def null(p): return G * np.log2(np.where(p % 4 == 1, 6.0, 2.0) * p) / p
def t1(p, L):
    return G * (np.log2(2 * p) + 0.5 * np.log2(3)) / p * np.interp(np.log(p), sizes, LF[:, LANES.index(L)])
for lo, hi, tag in ((5, B, 'all 50 from 5 to 136,279,841'), (1e6, B, 'above 1 million'), (1e7, B, 'above 10 million')):
    e = np.array([x for x in E if lo <= x < hi])
    q = P[(P >= lo) & (P < hi)]
    obs = [int(np.sum(e % 12 == L)) for L in LANES]
    nl = np.array([null(q[q % 12 == L]).sum() for L in LANES]); nl *= len(e) / nl.sum()
    tl = np.array([t1(q[q % 12 == L], L).sum() for L in LANES]); tl *= len(e) / tl.sum()
    cnt = np.array([np.sum(q % 12 == L) for L in LANES])
    print(f'{tag}: {len(e)} exponents')
    print('   class mod 12      ' + ''.join(f'{L:>9}' for L in LANES))
    print('   exponents found   ' + ''.join(f'{x:>9}' for x in obs))
    print('   plain heuristic   ' + ''.join(f'{x:>9.1f}' for x in nl))
    print('   theory 1 (lanes)  ' + ''.join(f'{x:>9.1f}' for x in tl))
    print('   all primes, share ' + ''.join(f'{x / cnt.sum():>9.5f}' for x in cnt) + f'   (class 1 behind the mean by {cnt.mean() - cnt[0]:,.0f} of {cnt.mean():,.0f})')
print('class of each exponent above 1 million: ' + ', '.join(f'{x:,}:{x % 12}' for x in E if x >= 1e6))
