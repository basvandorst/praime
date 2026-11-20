"""Rabbit clocks: for each prime p, z(p) = the least n >= 1 with p | F_n (rank of apparition; it divides p - (5/p) for
p != 5).  A prime divides some Lucas number L_n exactly when z(p) is even (L_n = F_2n / F_n).  Lagarias (1985, with an
erratum 1997): the primes dividing some Lucas number have density 2/3.  Also the share of p with z(p) = p - (5/p)
(maximal rank), an Artin-type question.  Prediction written before the run: the even share is within 0.005 of 2/3 at
10^6.   usage: rabbits.py X out.txt"""
import sys, math
from sympy import primerange, factorint, legendre_symbol
X = int(sys.argv[1]); out = open(sys.argv[2], 'w')
def pr(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')
def fib(n, m):                     # (F_n, F_{n+1}) mod m by fast doubling
    a, b = 0, 1
    for bit in bin(n)[2:]:
        c = a * ((2 * b - a) % m) % m; d = (a * a + b * b) % m
        a, b = (d, (c + d) % m) if bit == '1' else (c, d)
    return a
def z(p):
    if p == 2: return 3
    if p == 5: return 5
    n = p - legendre_symbol(5, p)
    for q, e in factorint(n).items():
        for _ in range(e):
            if fib(n // q, p) == 0: n //= q
            else: break
    return n
marks = [10 ** e for e in range(3, int(math.log10(X)) + 1)]
ev = 0; mx = 0; tot = 0; res = {}
for p in primerange(2, X + 1):
    zz = z(p); tot += 1; ev += (zz % 2 == 0)
    if p not in (2, 5) and zz == p - legendre_symbol(5, p): mx += 1
    for mk in marks:
        if p <= mk: res[mk] = (tot, ev, mx)
pr('X = %d: z(2) = %d, z(3) = %d, z(7) = %d, z(11) = %d (F_n: 3 | F_4, 7 | F_8, 11 | F_10)' % (X, z(2), z(3), z(7), z(11)))
pr('     x   primes   z(p) even (divides a Lucas number)   share   (2/3 = 0.6667)   z(p) = p - (5/p): share')
for mk in marks:
    t, e, m = res[mk]; pr('%8.0e %8d %20d %20.4f %22.4f' % (mk, t, e, e / t, m / t))
