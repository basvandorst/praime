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
ev = 0; mx = 0; tot = 0; res = {}; luc = 0; spl = 0; ine = 0; mod4 = 0
def lucas_div(p):                  # brute force: does p divide some L_n, n >= 1 (L_n mod p is periodic, period <= 2p + 2... use 6p + 6)
    a, b = 2 % p, 1 % p
    for n in range(1, 6 * p + 8):
        a, b = b, (a + b) % p
        if a == 0: return True
    return False
for p in primerange(2, X + 1):
    zz = z(p); tot += 1; ev += (zz % 2 == 0); luc += (zz % 2 == 0 or p == 2)
    if p < 3000: assert lucas_div(p) == (zz % 2 == 0 or p == 2), p
    if p not in (2, 5) and zz == p - legendre_symbol(5, p):
        mx += 1; mod4 += (p % 4 == 3)
        if legendre_symbol(5, p) == 1: spl += 1
        else: ine += 1
    for mk in marks:
        if p <= mk: res[mk] = (tot, ev, mx, luc, spl, ine, mod4)
pr('X = %d: z(2) = %d, z(3) = %d, z(7) = %d, z(11) = %d (F_n: 3 | F_4, 7 | F_8, 11 | F_10)' % (X, z(2), z(3), z(7), z(11)))
pr('checked for p < 3000 by brute force: p divides some L_n (n >= 1) exactly when z(p) is even or p = 2 (2 | L_3 = 4)')
pr('     x   primes   z(p) even   divide a Lucas number   share   (2/3 = 0.6667)   longest clock z(p) = p - (5/p) among p != 2, 5: count, share; split (z = p - 1), inert (z = p + 1); how many = 3 mod 4')
A = 0.3739558136
for mk in marks:
    t, e, m, l, sp, ie, m4 = res[mk]; pr('%8.0e %8d %10d %14d %20.4f %10d %10.4f %8d %8d %8d' % (mk, t, e, l, l / t, m, m / (t - 2), sp, ie, m4))
t, e, m, l, sp, ie, m4 = res[marks[-1]]
pr('at %d: split share %.5f (GRH value 9A/19 = %.5f), inert share %.5f (10A/19 = %.5f), total %.5f (Artin A = %.5f)' % (marks[-1], sp / (t - 2), 9 * A / 19, ie / (t - 2), 10 * A / 19, m / (t - 2), A))
