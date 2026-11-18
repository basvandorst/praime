"""Orders of the three-cell wheel (r = 7, u = 2 + y) modulo primes p that are not +-1 mod 7, and a search for numbers built from them."""
import sys, itertools, json
from sympy import primerange, factorint
from three_hunt_lib import mul, upow
from fold_check import fold_vec
P = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
pool = {}
for p in primerange(3, P):
    if p % 7 in (0, 1, 6): continue
    m = p*p + p + 1
    assert upow(m, p) == [1, 0, 0], p
    d = m
    for q in factorint(m):
        while d % q == 0 and upow(d//q, p) == [1, 0, 0]: d //= q
    pool[p] = d
small = sorted((d, p) for p, d in pool.items())
print('primes below', P, 'that are not +-1 mod 7:', len(pool), '| ten smallest periods (period, prime):', small[:10])
print('primes whose period is below the prime itself:', sum(1 for p, d in pool.items() if d < p), '| below 1000:', sum(1 for d in pool.values() if d < 1000))
cl = lambda x: {1: 0, 6: 0, 2: 1, 5: 1, 4: 2, 3: 2}[x % 7]
def ok_at(p, n):
    j = (cl(n) * (1 if cl(p) == 1 else 2)) % 3
    return j != 0 and (n - p**j) % pool[p] == 0
S = [p for d, p in small if d < 20000]
hits = []
for a, b in itertools.combinations(S, 2):
    n = a*b
    if n % 7 not in (1, 6) and ok_at(a, n) and ok_at(b, n): hits.append((n, (a, b)))
for a, b, c in itertools.combinations(S[:220], 3):
    n = a*b*c
    if n % 7 not in (1, 6) and ok_at(a, n) and ok_at(b, n) and ok_at(c, n): hits.append((n, (a, b, c)))
print('primes with a period below 20000:', len(S), '| products of two of them, or of three of the first 220, that meet the condition at every factor:', len(hits))
ex = []
for n, combo in sorted(hits)[:10]:
    okv = fold_vec(n, 7, 3)[0]; ex.append({'n': n, 'factors': combo, 'periods': [pool[p] for p in combo], 'home': okv})
    print(n, '=', ' * '.join(map(str, combo)), '| periods', [pool[p] for p in combo], '| three-cell wheel home, plain Python on the full ring:', okv)
json.dump({'P': P, 'pool': len(pool), 'small': len(S), 'hits': len(hits), 'examples': ex}, open('three_orders.json', 'w'), indent=1)
