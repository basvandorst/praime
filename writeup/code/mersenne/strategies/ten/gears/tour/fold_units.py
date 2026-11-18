"""Folding Pascal's ring of r cells (r prime) by a subgroup H of the units mod r: exact checks.
u_H = product over h in H of (x^h + 1) in Z[x]/(x^r - 1). H = the e-th powers mod r, e a prime power dividing (r - 1)/2.
1. Every u_H is constant on the classes of H (so it lives on e cells), and the product of all x^h + 1 is 1 mod Phi_r (a unit).
2. u_H is 1 mod Phi_r (the fold says nothing) exactly when 2 is an e-th power mod r.
3. r = 5: u_H = 2 + x + x^4 = phi^2, phi the golden ratio. So a number that comes home on 5 cells satisfies
   phi^(2n) = phi^2 (n = 1, 4 mod 5) or phi^(-2) (n = 2, 3 mod 5), hence n divides F(n - 1) or F(n + 1):
   it is a mistake of the Fibonacci wheel too. Checked on every liar of the 5-cell ring below 1e10."""
import json, cmath, math
from fold_check import prim_root, classes, folds_of
def polymul(a, b, r):
    c = [0]*r
    for i in range(r):
        if a[i]:
            for j in range(r):
                if b[j]: c[(i+j) % r] += a[i]*b[j]
    return c
def unit(r, e):
    cl = classes(r, e); P = [0]*r; P[0] = 1
    for h in range(1, r):
        if cl[h] == 0:
            f = [0]*r; f[0] = 1; f[h] += 1; P = polymul(P, f, r)
    assert all(P[z] == P[w] for z in range(1, r) for w in range(1, r) if cl[z] == cl[w])
    rep = {}
    for x in range(r-1, 0, -1): rep[cl[x]] = x
    return [P[rep[k]] - P[0] for k in range(e)], cl
primes = [p for p in range(5, 100) if all(p % d for d in range(2, p))]
bad = 0; rows = []; nothing = []
for r in primes:
    full = [0]*r; full[0] = 1
    for h in range(1, r):
        f = [0]*r; f[0] = 1; f[h] += 1; full = polymul(full, f, r)
    assert all(full[z] - full[0] == -1 for z in range(1, r)) or all(full[z] == full[1] for z in range(1, r))
    one = all(full[z] - full[0] == -1 for z in range(1, r))
    for e in folds_of(r):
        u, cl = unit(r, e)
        trivial = all(v == -1 for v in u)
        g = prim_root(r); ind = {}; x = 1
        for k in range(r-1): ind[x] = k; x = x*g % r
        two_in_H = ind[2] % e == 0
        if trivial != two_in_H: bad += 1
        if trivial: nothing.append((r, e))
        rows.append((r, e, trivial))
    if not one: bad += 1
print('primes r from 5 to 97:', len(primes), '| folds (r, e):', len(rows), '| product of all x^h + 1 is 1 mod Phi_r for every r:', bad == 0 or 'see below')
print('folds that say nothing (u = 1):', nothing)
print('"says nothing" happens exactly when 2 is an e-th power mod r; exceptions:', bad)
u5, cl5 = unit(5, 2)
z = cmath.exp(2j*cmath.pi/5); val = ((1+z)*(1+z**4)).real; phi = (1+5**0.5)/2
print('r = 5, e = 2: u in the basis (x + x^4, x^2 + x^3) =', u5, '| as a number', round(val, 12), '| phi^2 =', round(phi*phi, 12))
def fib_pair(n, m):
    if n == 0: return (0, 1)
    a, b = fib_pair(n >> 1, m)
    c = a*((2*b - a) % m) % m; d = (a*a + b*b) % m
    return (d, (c + d) % m) if n & 1 else (c, d)
L5 = [n for n in json.load(open('ring_psp2_to1e10.json'))['5'] if n % 2]
ok = sum(1 for n in L5 if fib_pair(n - 1 if n % 5 in (1, 4) else n + 1, n)[0] == 0)
print('liars of the 5-cell ring below 1e10:', len(L5), '| of these also liars of the Fibonacci wheel (n divides F(n -+ 1)):', ok)
def phi2pow(n):
    a, b = 1, 0
    def mul(p, q): return ((p[0]*q[0] + p[1]*q[1]) % n, (p[0]*q[1] + p[1]*q[0] + p[1]*q[1]) % n)
    base = (1 % n, 1 % n); res = (1 % n, 0); e = n
    while e:
        if e & 1: res = mul(res, base)
        base = mul(base, base); e >>= 1
    return res
cnt = imp = 0; comp = []
from sympy import isprime
for n in range(7, 200000, 2):
    if n % 5 in (2, 3):
        if phi2pow(n) == (2 % n, n - 1):
            cnt += 1
            if fib_pair(n + 1, n)[0] == 0: imp += 1
            if not isprime(n): comp.append(n)
np_ = sum(1 for p in range(7, 200000, 2) if p % 5 in (2, 3) and isprime(p))
print('odd n below 200000 with n = 2, 3 mod 5 and (phi^2)^n = phi^(-2) mod n:', cnt, '| of these n divides F(n + 1):', imp, '| primes among them:', cnt - len(comp), 'of', np_, 'primes in those classes | not prime:', len(comp), 'first', comp[:4])
u7 = 2 + 2*math.cos(2*math.pi/7)
print('r = 7, e = 3: u = 2 + 2 cos(2 pi / 7) =', round(u7, 10), '| u^3 - 5u^2 + 6u - 1 =', round(u7**3 - 5*u7**2 + 6*u7 - 1, 9))
