"""Round 18 check: n^2+n+41. Counts of prime values in F1, F2, F3; BH constant; pair series for d = 1 (collision rule)."""
import math, json
import numpy as np, sympy
from sympy.ntheory import legendre_symbol, sqrt_mod
f = lambda n: n * n + n + 41
def roots(p):
    if p == 2: return []
    return sorted(set(n for n in sqrt_mod(-163 % p, p, all_roots=True) for n in [(n - 1) * pow(2, -1, p) % p]))
# verify the root rule on small primes by brute force
for p in list(sympy.primerange(2, 400)):
    br = [n for n in range(p) if f(n) % p == 0]
    assert br == roots(p), p
    if p < 41: assert br == []
print('roots rule OK for p < 400; no roots for p < 41')
def count(L, W, B=200000):
    M = np.ones(W, dtype=bool)
    for p in sympy.primerange(3, B):
        for r in roots(p):
            st = (r - L) % p
            M[st::p] = False
    c = 0
    for i in np.nonzero(M)[0]:
        if sympy.isprime(f(L + int(i))): c += 1
    return c
for name, L, W in (('F1', 10 ** 6, 10 ** 6), ('F2', 10 ** 8, 5 * 10 ** 5), ('F3', 10 ** 10, 10 ** 5)):
    print(name, 'prime values', count(L, W), flush=True)
# BH constant and pair series (d = 1), product to P
PR = list(sympy.primerange(2, 20_000_000))
def chi(p): return 0 if p == 163 else legendre_symbol(-163 % p, p)
lc = 0.0; lp = 0.0; lp2000 = None
for p in PR:
    w = 0 if p == 2 else 1 + chi(p)
    lc += math.log(1 - w / p) - math.log(1 - 1 / p)
    # pair f(n) f(n+1): |R u (R-1)|
    if p == 2: w2 = 0
    elif p == 163: w2 = 2
    else: w2 = 2 * w - (1 if (w == 2 and (1 + 163) % p == 0) else 0)
    lp += math.log(1 - w2 / p) - 2 * math.log(1 - 1 / p)
    if p < 2000: lp2000 = lp
    if p in (100003, 1000003):
        print('partial at', p, 'C =', math.exp(lc), 'pair =', math.exp(lp))
print('BH constant to 2e7: %.5f (literature 2 x 3.3197732 = 6.63955)' % math.exp(lc))
print('pair series d=1 to 2e7: %.4f; to 2000: %.4f' % (math.exp(lp), math.exp(lp2000)))
# brute-force check of omega_2 at small p
for p in list(sympy.primerange(3, 300)):
    br = len([n for n in range(p) if f(n) * f(n + 1) % p == 0])
    w = 1 + chi(p) if p != 163 else 1
    w2 = 2 * w - (1 if (w == 2 and 164 % p == 0) else 0)
    assert br == w2, p
print('collision rule OK for d = 1, p < 300')
