"""Why the digit-sum rules lean: the same lean in all numbers without the factor 3 (binary) or 11 (decimal)."""
import math
from sympy import primerange
def ds(n, b):
    s = 0
    while n: s += n % b; n //= b
    return s
X = 10 ** 6
P = list(primerange(2, X + 1))
for b, cases in ((2, ((2,), (2, 3))), (10, ((2, 3, 5), (2, 3, 5, 11), (2, 3, 5, 7, 13), (2, 3, 5, 7, 11, 13)))):
    for mods in cases:
        for top in (5000, X):
            o = e = 0
            for n in range(1, top + 1):
                if any(n % m == 0 for m in mods): continue
                if ds(n, b) % 2: o += 1
                else: e += 1
            print('base %d, numbers up to %d without the factors %s: odd digit sum %.4f' % (b, top, mods, o / (o + e)))
    for top in (5000, X):
        Q = [p for p in P if p <= top]
        o = sum(ds(p, b) % 2 for p in Q)
        print('base %d, the %d primes up to %d: odd digit sum %.4f (one standard deviation of chance: %.4f)' % (b, len(Q), top, o / len(Q), 0.5 / math.sqrt(len(Q))))
