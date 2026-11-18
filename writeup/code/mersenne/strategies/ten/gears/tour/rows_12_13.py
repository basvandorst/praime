"""The three results of machines 12 and 13 for 250 wheels (lamps stay on), in plain Python:
the three-cell wheel on the full ring (fold_check.fold_vec) and Wilson's wheel ((n - 1)! = -1 mod n)."""
import json
from sympy import isprime
from fold_check import fold_vec, classes
SIZES = [7, 13, 19, 37, 61, 67, 73, 79, 97, 103, 139, 151, 163, 181, 193, 199, 211, 241]
def three(n):
    if n == 2: return True
    if n % 2 == 0: return False
    for r in SIZES:
        if n % r == 0:
            if n == r: continue
            return False
        if classes(r, 3)[n % r] != 0: return fold_vec(n, r, 3)[0]
    return False
def wilson(n):
    f = 1
    for k in range(2, n): f = f*k % n
    return f == n - 1
MX = {2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127}; out = {}
for name, home in (('three', three), ('wilson', wilson)):
    lit = []; best = {'p': (0, 0), 'm': (0, 0), 'w': (0, 0)}
    for t in range(2, 251):
        if home(t): lit.append(t)
        p = [k for k in lit if isprime(k)]; m = [k for k in p if k in MX]; w = [k for k in p if isprime(k - 2) or isprime(k + 2)]
        for key, a in (('p', p), ('m', m), ('w', w)):
            if len(a) > best[key][0]: best[key] = (len(a), t)
    row = list(best['p'] + best['m'] + best['w']); out[name] = row
    print("%s, 250 wheels: primes %d (tick %d), Mersenne exponents %d (tick %d), twin primes %d (tick %d); lamps of numbers that are not prime: %d" % ((name,) + tuple(row) + (sum(1 for k in lit if not isprime(k)),)))
    assert all(isprime(k) for k in lit)
json.dump(out, open('rows_12_13.json', 'w'))
