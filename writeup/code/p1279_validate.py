"""Puzzle 1279 Q3: least Sierpinski and Riesel numbers whose covering sets need 12 primes.

Independent checks in Python (plain integers, numpy for the list check):
  1. data/p1279_primes.txt is exactly the set of primes p < 2^32 with ord_p(2) | 10080
     (small orders: full factorisation of 2^d - 1; large orders: every p = 1 + j*d < 2^32 tested
     with a vectorised 2^d mod p);
  2. for each claimed number: the claimed set covers k*2^n + s for every n, it is minimal, and the
     minimal covers drawn from the list are exactly the claimed ones (all of size 12);
  3. the same for the known terms a(6) .. a(11) of A206001 / A206430 (sizes 6 .. 11);
  4. the numbers with three or more minimal covers listed in A394518 / A394519 below the search
     bounds have three or more covers here too.
Run: python3 p1279_validate.py   (about 2 minutes)
"""
import os
from math import gcd
import numpy as np
from sympy import factorint, isprime, divisors

L = 10080
here = os.path.dirname(os.path.abspath(__file__))
lam = [tuple(map(int, l.split())) for l in open(os.path.join(here, '..', 'data', 'p1279_primes.txt'))]

def order2(p):
    t = p - 1
    for q in factorint(p - 1):
        while t % q == 0 and pow(2, t // q, p) == 1: t //= q
    return t

# 1. the prime list
assert all(isprime(p) and order2(p) == m and L % m == 0 for p, m in lam)
want = set()
for d in divisors(L):
    if 2 <= d <= 64:
        want |= {p for p in factorint(2**d - 1) if p < 2**32 and order2(p) == d}
    elif d > 64:
        step = d if d % 2 == 0 else 2 * d
        for lo in range(1, 2**32, step * 2**22):
            p = np.arange(lo + step, min(lo + step * (2**22 + 1), 2**32), step, dtype=np.uint64)
            r = np.ones_like(p); b = np.full_like(p, 2); e = d
            while e:                                    # 2^d mod p, all products < 2^64
                if e & 1: r = (r * b) % p
                b = (b * b) % p; e >>= 1
            for x in p[r == 1].tolist():
                if isprime(x) and order2(x) == d: want.add(x)
assert want == {p for p, m in lam}, (sorted(want - {p for p, m in lam}), sorted({p for p, m in lam} - want))
print(f"prime list: {len(lam)} primes p < 2^32 with ord_p(2) | {L}, complete")

FULL = (1 << L) - 1
def cls_mask(a, m):
    x = 0
    for n in range(a, L, m): x |= 1 << n
    return x

def minimal_covers(k, s):
    used = []
    for p, m in lam:
        a = [a for a in range(m) if (k * pow(2, a, p) + s) % p == 0]
        if a: used.append((p, cls_mask(a[0], m)))
    covers = set()
    def rec(i, chosen, cov):
        if cov == FULL:
            for j in chosen:       # minimal: dropping any prime uncovers a point
                o = 0
                for t in chosen:
                    if t != j: o |= used[t][1]
                if o == FULL: return
            covers.add(tuple(sorted(used[j][0] for j in chosen)))
            return
        if len(chosen) >= 20 or i == len(used): return
        # include used[i] only if it adds a point
        if used[i][1] & ~cov: rec(i + 1, chosen + [i], cov | used[i][1])
        rest = cov
        for t in range(i + 1, len(used)): rest |= used[t][1]
        if rest == FULL: rec(i + 1, chosen, cov)
    rec(0, [], 0)
    return sorted(covers, key=lambda c: (len(c), c))

def covers_all_n(k, s, cs):
    M = 1
    for p in cs:
        m = order2(p); M = M * m // gcd(M, m)
    return all(any((k * pow(2, n, p) + s) % p == 0 for p in cs) for n in range(M))

# 2. the new numbers
claims = [(228288720707, 1, (3, 5, 7, 11, 13, 17, 31, 37, 41, 151, 181, 433)),
          (22828394579, -1, (3, 5, 7, 11, 13, 17, 31, 37, 151, 331, 433, 631))]
for k, s, cs in claims:
    assert covers_all_n(k, s, cs) and k > max(cs)
    for p in cs: assert not covers_all_n(k, s, tuple(q for q in cs if q != p))
    mc = minimal_covers(k, s)
    assert mc == [cs], mc
    print(f"{'Sierpinski' if s > 0 else 'Riesel'} {k}: unique minimal cover {{{', '.join(map(str, cs))}}} (12 primes)")

# 3. known terms
known = {1: [271129, 78557, 1777613, 169073869, 11679013547, 1143502909],
         -1: [509203, 777149, 7106977, 60014203, 10909924841, 10591451383]}
for s, ks in known.items():
    for n, k in enumerate(ks, start=6):
        sizes = sorted({len(c) for c in minimal_covers(k, s)})
        assert n in sizes, (k, sizes)
    print(f"A{'206001' if s > 0 else '206430'} a(6..11): sizes confirmed")

# 4. numbers with at least three covering sets (OEIS A394518, A394519) inside the search ranges
for s, ks, bound in ((1, [56889499253, 105839880557, 111056109533, 129119577739], 228288720707),
                     (-1, [14995746487], 22828394579)):
    for k in ks:
        assert k < bound and len(minimal_covers(k, s)) >= 3
print("A394518 / A394519 terms below the bounds: three or more minimal covers each")

# 5. every number with three or more minimal covers found by the search (data/p1279_multi.txt):
#    each listed set covers and is minimal, and the list equals the minimal covers from the prime list
n5 = 0; multi = {1: {}, -1: {}}
for line in open(os.path.join(here, '..', 'data', 'p1279_multi.txt')):
    if line.startswith('#') or not line.strip(): continue
    f = line.strip().split(' | ')
    side, k, cnt = f[0].split(); s = 1 if side == 'S' else -1; k = int(k); cnt = int(cnt)
    sets = sorted((tuple(sorted(map(int, x.split()))) for x in f[1:]), key=lambda c: (len(c), c))
    assert len(sets) == cnt
    for cs in sets:
        assert covers_all_n(k, s, cs)
        for p in cs: assert not covers_all_n(k, s, tuple(q for q in cs if q != p))
    assert minimal_covers(k, s) == sets, k
    n5 += 1; multi[s][k] = sets
print(f"{n5} numbers with three or more minimal covers: every set checked, and the lists are complete")

# 6. the least numbers with 3, 4, 5, 6 covering sets (the data file holds every k with >= 3 below the bounds)
for s, name in ((1, 'Sierpinski'), (-1, 'Riesel')):
    M = multi[s]
    for c in (3, 4, 5, 6):
        ge = [k for k in M if len(M[k]) >= c]; eq = [k for k in M if len(M[k]) == c]
        print(f"   {name}: least with >= {c} sets: {min(ge) if ge else 'none'}"
              f" ({len(M[min(ge)]) if ge else '-'} sets); least with exactly {c}: {min(eq) if eq else 'none'}")
assert min(multi[1]) == 3278650021 and len(multi[1][3278650021]) == 6
assert not [k for k in multi[1] if len(multi[1][k]) == 5]
assert min(k for k in multi[-1] if len(multi[-1][k]) >= 5) == 141719191807
assert min(k for k in multi[-1] if len(multi[-1][k]) == 5) == 203675118523
# the OEIS lists A394518 / A394519 are exactly the k with three sets of at most 7 primes that avoid 673
oeis = {1: [56889499253, 105839880557, 111056109533, 129119577739],
        -1: [14995746487, 28299882787, 152078237473, 158736554201, 171366858989, 185819859341]}
for s in (1, -1):
    r = sorted(k for k, sets in multi[s].items() if sum(1 for c in sets if len(c) <= 7 and 673 not in c) >= 3)
    assert r == oeis[s], r
print("A394518 / A394519 below the bounds = the k with three covering sets of <= 7 primes avoiding 673")
print("all checks passed")
