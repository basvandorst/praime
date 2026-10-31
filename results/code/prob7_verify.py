"""Problem 7: parity of x and its inverse mod p. Checks the counts for x = 2..p-1 (x = 1 left out, as on the page)
and the involution x -> p - x behind Np(oo) = Np(ee) - 1, for all odd primes p < 5000 and odd composites n < 2000."""
from sympy import primerange, isprime
from math import gcd
def counts(n):
    c = {'ee': 0, 'oo': 0, 'op': 0}
    for x in range(2, n):
        if gcd(x, n) != 1: continue
        y = pow(x, -1, n)
        assert (n - x) * (n - y) % n == 1                     # (-x)(-x') = xx' = 1
        assert (n - x) % 2 != x % 2 and (n - y) % 2 != y % 2  # n odd: n - x flips parity
        c['ee' if x % 2 == 0 and y % 2 == 0 else 'oo' if x % 2 and y % 2 else 'op'] += 1
    return c
for p in primerange(3, 5000):
    c = counts(p)
    assert sum(c.values()) == p - 2 and c['oo'] == c['ee'] - 1 and c['op'] % 2 == 0, p
for n in range(9, 2000, 2):
    if isprime(n): continue
    c = counts(n)
    assert c['oo'] == c['ee'] - 1, n                          # the same holds for every odd n
assert counts(13)['op'] == 6
# Lehmer's question: Np(op) - (p-1)/2 is of size about sqrt(p) (Zhang Wenpeng 1994: O(p^(1/2) log^2 p))
worst = max(abs(counts(p)['op'] - (p - 1) / 2) / p**0.5 for p in primerange(1000, 5000))
print(f"max |Np(op) - (p-1)/2| / sqrt(p) for 1000 < p < 5000: {worst:.2f}")
print("all checks passed")
