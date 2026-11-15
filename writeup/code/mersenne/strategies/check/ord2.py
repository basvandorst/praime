"""Independent check of the ord2 note with sympy.n_order: odd primes q <= 500,000 whose order of 2 is prime, index 1 share,
the even-index law for q = +-1 mod 8, and the indices of 8,191 and 131,071."""
import sympy
tot = pf = i1 = odd = 0
for q in sympy.primerange(3, 500001):
    d = sympy.n_order(2, q); i = (q - 1) // d; tot += 1
    if sympy.isprime(d): pf += 1
    if i == 1: i1 += 1
    if i % 2 == 1 and q % 8 in (1, 7): odd += 1
print(f'odd primes {tot}; prime order of 2: {pf} ({100*pf/tot:.2f}%); index 1: {i1} ({100*i1/tot:.2f}%); odd index with q = +-1 mod 8: {odd}')
print('index of 8191:', 8190 // sympy.n_order(2, 8191), '; of 131071:', 131070 // sympy.n_order(2, 131071))
from collections import Counter
c = Counter(); bad = 0
for q in sympy.primerange(5, 500001):
    d = sympy.n_order(2, q)
    if sympy.isprime(d) and d > 2:
        i = (q - 1) // d; c[i % 8] += 1
        if i % 2 or i % 8 == 4: bad += 1
print('Mersenne factors q (prime order d > 2) by index mod 8:', dict(sorted(c.items())), '; odd or 4 mod 8:', bad)
