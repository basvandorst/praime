"""Referee check S5: local means of the class-number factor F_l(D), D = 4pN - r^2 N^2, over p in Z_l^x.
(1) exact enumeration of p mod l^J with an enclosure for undetermined residues;
(2) direct check over actual primes p (coredisc): for every l | N the conductor f0 is prime to l and chi(l) = 0.
F_l(D) := sum_{j<a} l^-j + l^-a/(1 - chi(l)/l), -D = Delta0 f0^2, a = v_l(f0)  (local factor of pi H(D)/sqrt(D))."""
import sys
from fractions import Fraction as Fr
import sympy as sp
import cypari2
pari = cypari2.Pari()

def Floc(l, a, chi):
    return sum(Fr(1, l ** j) for j in range(a)) + Fr(1, l ** a) / (1 - Fr(chi, l))

def local_from_residue(l, J, m):
    """m = -D mod l^J (as integer in [0, l^J)). Returns F or None if undetermined."""
    M = l ** J
    m %= M
    if m == 0:
        return None
    v = 0
    while m % l == 0:
        m //= l; v += 1
    # unit part u = m mod l^(J-v)
    if l != 2:
        if v + 1 > J:
            return None
        a = v // 2
        chi = 0 if v % 2 else sp.legendre_symbol(m % l, l)
        return Floc(l, a, chi)
    if v + 3 > J:
        return None
    u8 = m % 8
    if v % 2 == 0:
        b = v // 2
        if u8 % 4 == 1:
            return Floc(2, b, 1 if u8 == 1 else -1)
        if b == 0:
            raise ValueError('not a discriminant')
        return Floc(2, b - 1, 0)
    b = (v - 1) // 2
    if b == 0:
        raise ValueError('not a discriminant')
    return Floc(2, b - 1, 0)

def e_gen(l, r):
    u = Fr(1, (l - 1) ** 2 * (l + 1))
    return Fr(l * l, l * l - 1) if r % l == 0 else 1 - u

Js = {2: 13, 3: 8, 5: 6, 7: 5}
levels = [20014, 20019, 20058, 20015, 10006, 10014, 20021, 20023, 20027, 7 * 3343, 2 * 3 * 5 * 7 * 11 * 13 * 17]
levels = [N for N in levels if sp.factorint(N) and all(e == 1 for e in sp.factorint(N).values())]
maxdev = 0
ncase = 0
for N in levels:
    for l, J in Js.items():
        M = l ** J
        for r in range(0, 9):
            tot = Fr(0); und = 0; cnt = 0
            for p in range(1, M):
                if p % l == 0:
                    continue
                cnt += 1
                F = local_from_residue(l, J, -(4 * p * N - r * r * N * N))
                if F is None:
                    und += 1
                else:
                    tot += F
            lo = float(tot / cnt); hi = float((tot + und * Fr(l, l - 1)) / cnt)
            pred = 1 if N % l == 0 else e_gen(l, r)
            ok = lo - 1e-12 <= float(pred) <= hi + 1e-12
            ncase += 1
            maxdev = max(maxdev, hi - lo)
            if not ok:
                print('MISMATCH', N, l, r, lo, hi, float(pred))
print(f'enumeration: {ncase} cases (levels {levels}), all inside enclosure unless printed; max width {maxdev:.2e}')

# (2) actual primes: for l | N, v_l(f0) = 0 and kron(Delta0, l) = 0 for all p <= 4N, r < 2 sqrt(p/N)
for N in [20014, 20058, 10006, 20019, 20015]:
    bad = 0; tot = 0
    qs = list(sp.factorint(N))
    for p in sp.primerange(3, 4 * N):
        if N % p == 0:
            continue
        r = 0
        while r * r * N < 4 * p:
            D = 4 * p * N - r * r * N * N
            d0 = int(pari.coredisc(-D))
            f2 = (-D) // d0
            for q in qs:
                if f2 % q == 0 or int(pari.kronecker(d0, q)) != 0:
                    bad += 1
            # 2-adic valuation of D for N even (diagnostic)
            tot += 1
            r += 1
    print(f'N={N} factors {qs}: {tot} discriminants, violations of (l | N => l !| f0, chi(l)=0): {bad}')
# valuation of D at 2 for N even: distribution
N = 20014
from collections import Counter
c = Counter()
for p in sp.primerange(3, 4 * N):
    for r in range(0, 4):
        if r * r * N < 4 * p:
            D = 4 * p * N - r * r * N * N
            v = (D & -D).bit_length() - 1
            c[(r % 2, v)] += 1
print('N = 20014: (r mod 2, v_2(D)) counts:', dict(c))
