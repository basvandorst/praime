"""Independent checks (plain Python pow, sympy) of the wieferich and cunningham notes, N = 500,000."""
import sympy, math
N = 500000
P = list(sympy.primerange(3, N + 1))
qs = {p: ((pow(2, p - 1, p * p) - 1) // p) / p for p in P}
wf = [p for p in P if qs[p] == 0]
ys = sorted(qs.values()); n = len(ys)
ks = max(max(abs(y - i / n), abs(y - (i + 1) / n)) for i, y in enumerate(ys))
near = min((v, p) for p, v in qs.items() if v > 0)
print('Wieferich primes:', wf, '; orders of 2:', [sympy.n_order(2, p) for p in wf], [sympy.factorint(sympy.n_order(2, p)) for p in wf])
print(f'odd primes {n}; top half {100*sum(1 for y in ys if y < .5)/n:.2f}%; KS {ks:.4f} (5% critical 1.36/sqrt(n) = {1.36/math.sqrt(n):.4f}); nearest {near[1]} at {near[0]:.2e}')
def L(p):
    l = 1
    while sympy.isprime(2 * p + 1): p = 2 * p + 1; l += 1
    return l
lens = {p: L(p) for p in sympy.primerange(2, N + 1)}
from collections import Counter
c = Counter(lens.values()); sg = sum(1 for v in lens.values() if v > 1)
pred = 1.3203236316937391 * sum(1 / (math.log(k) * math.log(2 * k + 1)) for k in range(3, N + 1))
E = [2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243,110503,132049,216091]
print('chain lengths:', dict(sorted(c.items())), '; Sophie Germain:', sg, '; HL prediction', round(pred), '; length-6 starts:', [p for p, v in lens.items() if v == 6])
print('exponents with 2p+1 prime:', [p for p in E if lens[p] > 1], '; exponents = 3 mod 4 above 3 with 2p+1 prime:', [p for p in E if p > 3 and p % 4 == 3 and lens[p] > 1])
