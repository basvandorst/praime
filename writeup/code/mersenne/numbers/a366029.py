"""A366029: k > 1 with (k-1)^(2^k - 2) = 1 mod (2^k - 1) k^2. Direct check of the definition for k <= X,
against the split: (A) 2^k = 2 mod k [from the k^2 part] and (B) (k-1)^(2^k-2) = 1 mod 2^k - 1 [Fermat test of M_k, base k-1]."""
import sys, gmpy2
X = int(sys.argv[1])
E = {2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941}
direct, split = [], []
for k in range(2, X + 1):
    M = (gmpy2.mpz(1) << k) - 1
    if pow(gmpy2.mpz(k - 1), M - 1, M * k * k) == 1: direct.append(k)
    if pow(2, k, k) == 2 % k and pow(gmpy2.mpz(k - 1), M - 1, M) == 1: split.append(k)
print('direct definition, k <=', X, ':', direct)
print('split (A and B) equals direct:', direct == split)
print('extras (not Mersenne exponents):', [k for k in direct if k not in E], ' missing exponents:', [p for p in sorted(E) if p <= X and p not in direct])
print('extras of the form 2^m + 1:', [k for k in direct if k not in E and (k - 1) & (k - 2) == 0])
