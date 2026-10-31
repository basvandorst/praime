# Puzzle 819: every g with g(g(n)) = p(n) (the n-th prime) comes from a perfect matching of Pebody's rows plus an
# orientation per pair; we check f (Pebody), g (Vantieghem) and a random third solution on n <= NMAX.
import random
from sympy import prime, primepi, isprime
NMAX = 3000
def AB(n):                               # n = Z(a, b): row a, position b
    b = 1
    while isprime(n): n = int(primepi(n)); b += 1
    return n - int(primepi(n)), b
def Z1(m):                               # m-th non-prime (1 counts as non-prime)
    a = m + int(primepi(m))
    while a - primepi(a) != m: a += 1
    return a
def Z(a, b):
    x = Z1(a)
    for _ in range(b - 1): x = prime(x)
    return x
def make(partner, first):                # partner[a] = b; first[a] = True if g(Z(a,k)) = Z(b,k) (a is the "leading" row)
    def g(n):
        a, b = AB(n); c = partner[a]
        return Z(c, b) if first[a] else Z(c, b + 1)
    return g
R = 10000                                 # rows 1 .. 2R paired (enough for the values reached)
# f: pairs (2m+1, 2m+2), odd row leads;  Vantieghem's g: same pairs, even row leads
pairs_f = {a: (a + 1 if a % 2 else a - 1) for a in range(1, 2 * R + 1)}
f = make(pairs_f, {a: a % 2 == 1 for a in pairs_f})
g = make(pairs_f, {a: a % 2 == 0 for a in pairs_f})
random.seed(1); rows = list(range(1, 401)); random.shuffle(rows)          # random matching of the first 400 rows
pairs_h, first_h = dict(pairs_f), {a: a % 2 == 1 for a in pairs_f}       # standard pairing beyond
for i in range(0, len(rows), 2):
    x, y = rows[i], rows[i + 1]; pairs_h[x], pairs_h[y] = y, x
    lead = random.random() < 0.5; first_h[x], first_h[y] = lead, not lead
h = make(pairs_h, first_h)
for fun, name in ((f, "f"), (g, "g"), (h, "h")):
    assert all(fun(fun(n)) == prime(n) for n in range(1, 400)), name
assert [f(n) for n in range(1, 21)] == [4, 7, 17, 2, 59, 8, 3, 13, 10, 23, 277, 14, 19, 37, 16, 47, 5, 20, 41, 61]
assert f(2016) == 2018 and g(2016) == 17551
assert any(h(n) != f(n) and h(n) != g(n) for n in range(1, 50))
print("f, g and a random third solution all satisfy g(g(n)) = p(n) for n < 400")
print("all claims verified")
