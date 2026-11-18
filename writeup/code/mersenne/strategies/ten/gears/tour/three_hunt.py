"""Try to build a number that fools the three-cell wheel when it stands alone (no doubling wheel).
For r = 7 the wheel turns by u = 2 + y, y = 2 cos(2 pi / 7), y^3 + y^2 - 2y - 1 = 0. Let p be a prime that is not +-1 mod 7 (the three
cells stay one field mod p) and d the period of u mod p. A number n with the prime factor p is home mod p exactly when d divides
n - p^j, where j = 1 or 2 is fixed by n = +-p^j (mod 7). For most primes d is larger than p (three_orders.py), which rules out small
products. Primes with a short period d are the prime factors of Norm(u^d - 1), and they are huge. This script factors those norms with
PARI, collects the primes, and looks for products of two or three of them that meet the condition at every factor.
Every hit is checked on the full ring in plain Python (fold_check.fold_vec)."""
import sys, itertools, json
from math import gcd
import cypari2
from three_hunt_lib import mul, norm, upow, U
from fold_check import fold_vec
pari = cypari2.Pari(); pari.allocatemem(600*10**6)
D = int(sys.argv[1]) if len(sys.argv) > 1 else 80
pool = {}; A = {}; cur = [1, 0, 0]; skipped = []
for d in range(1, D + 1):
    cur = mul(cur, U); v = cur[:]; v[0] -= 1; N = abs(norm(v)); A[d] = N
    if any(q % 3 == 2 for q in range(2, d + 1) if d % q == 0 and all(q % t for t in range(2, q))): continue
    prim = N
    for k in range(1, d):
        if d % k == 0:
            g = gcd(prim, A[k])
            while g > 1: prim //= g; g = gcd(prim, A[k])
    if prim < 2: continue
    if prim > 10**66: skipped.append(d); continue
    F = pari.factor(prim)
    for i in range(len(F[0])):
        p = int(F[0][i])
        if p % 7 in (2, 3, 4, 5) and p != 2:
            dd = next(k for k in range(1, d + 1) if d % k == 0 and upow(k, p) == [1, 0, 0])
            if dd == d: pool[p] = d
ps = sorted(pool)
print('periods d up to', D, '(norms above 1e66 left out: d =', skipped, ') | primes not +-1 mod 7 with such a period:', len(ps))
print('  (prime digits, period):', [(len(str(p)), pool[p]) for p in ps])
cl = lambda x: {1: 0, 6: 0, 2: 1, 5: 1, 4: 2, 3: 2}[x % 7]
def ok_at(p, n):
    j = (cl(n) * (1 if cl(p) == 1 else 2)) % 3
    return j != 0 and (n - p**j) % pool[p] == 0
hits = []; tried = 0
for k in (2, 3, 4):
    for combo in itertools.combinations(ps, k):
        n = 1
        for p in combo: n *= p
        if n % 7 in (1, 6): continue
        tried += 1
        if all(ok_at(p, n) for p in combo): hits.append((n, combo))
print('products of 2, 3 or 4 of these primes that do not fit the ring of 7 cells:', tried, '| meet the condition at every factor:', len(hits))
ex = []
for n, combo in sorted(hits)[:8]:
    okv = fold_vec(n, 7, 3)[0]; ex.append({'n': str(n), 'factors': [str(p) for p in combo], 'periods': [pool[p] for p in combo], 'home': okv})
    print('n =', n, '(%d digits) =' % len(str(n)), ' * '.join(map(str, combo)), '| periods', [pool[p] for p in combo], '| three-cell wheel home, plain Python on the full ring:', okv, '| 2^n = 2 (mod n):', pow(2, n, n) == 2)
json.dump({'D': D, 'pool': [[str(p), pool[p]] for p in ps], 'tried': tried, 'hits': len(hits), 'examples': ex}, open('three_hunt.json', 'w'), indent=1)
