"""Independent check of round 103 (sympy, plain integers): windows of r consecutive prime values of n^2 + n + 41
for n in [1e3, 2e5), compared with the gmpy2 counts restricted to the same range; the longest run found by runs_data.py;
S_1 against the known Bateman-Horn constant of Euler's polynomial. usage: runs_check.py"""
import json, os
from sympy import isprime
import gmpy2

f = lambda n: n * n + n + 41
A, B = 10 ** 3, 2 * 10 ** 5
s = [isprime(f(n)) for n in range(A, B + 10)]
g = [bool(gmpy2.is_prime(f(n))) for n in range(A, B + 10)]
assert s == g
for r in range(1, 11):
    c = sum(all(s[i + j] for j in range(r)) for i in range(B - A))
    print(r, c)
d = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'runs_data.json')))
n0, L = d['longest_run']['start_n'], d['longest_run']['length']
assert all(isprime(f(n)) for n in range(n0, n0 + L)) and not isprime(f(n0 - 1)) and not isprime(f(n0 + L))
law = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'runs_law.json')))
assert abs(law['1']['S_r'] - 2 * 3.3197732) < 2e-3          # Bateman-Horn constant of n^2 + n + 41 (literature 6.63955)
# windows from the run lengths must equal the window counts
for r in range(1, 9):
    w = sum((int(k) - r + 1) * v for k, v in d['runs_by_length'].items() if int(k) >= r)
    assert w == d['rows'][str(r)]['observed'], r
print('longest run', L, 'from n =', n0, 'checked; window counts consistent with run lengths; all checks pass')
