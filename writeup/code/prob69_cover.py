# Problem 69: a Brier number whose + and - covering sets are disjoint.
import itertools, random, sys
from sympy import isprime, divisors, n_order
# + side: classic covering of n mod 24 with primes 3,5,7,13,17,241 (orders 2,4,3,12,8,24)
plus_primes = {3: 2, 5: 4, 7: 3, 13: 12, 17: 8, 241: 24}
def find_plus():
    ps = list(plus_primes)
    for res in itertools.product(*[range(plus_primes[p]) for p in ps]):
        if all(any(n % plus_primes[p] == r for p, r in zip(ps, res)) for n in range(24)):
            return dict(zip(ps, res))
plus = find_plus(); print("+ covering (prime: residue of n mod ord):", plus)
# - side: primes not in plus, with order m dividing L
L = int(sys.argv[1]) if len(sys.argv) > 1 else 5040
avail = {}
for m in divisors(L):
    if m < 5: continue
    found = []
    for k in range(1, 200000):
        p = k * m + 1
        if p > 10**8: break
        if p in plus_primes or not isprime(p): continue
        if pow(2, m, p) == 1 and n_order(2, p) == m: found.append(p)
        if len(found) >= 6: break
    if found: avail[m] = found
dens = sum(len(v) / m for m, v in avail.items())
print(f"L={L}: moduli with primes {len(avail)}, density sum {dens:.3f}")
# greedy covering of Z/L with restarts
best = None
random.seed(1)
items = [(m, p) for m, ps in avail.items() for p in ps]
for trial in range(300):
    covered = bytearray(L); choice = {}
    random.shuffle(items)
    left = L
    pool = sorted(items, key=lambda x: x[0])
    usedp = set()
    while left:
        bestgain, bm = 0, None
        for m, p in pool:
            if p in usedp: continue
            # best residue for this modulus
            for r in range(m):
                g = sum(1 for n in range(r, L, m) if not covered[n])
                if g > bestgain or (g == bestgain and bm and m < bm[0]): bestgain, bm = g, (m, p, r)
        if not bm or bestgain == 0: break
        m, p, r = bm; usedp.add(p); choice[p] = (m, r)
        for n in range(r, L, m):
            if not covered[n]: covered[n] = 1; left -= 1
    if left == 0:
        best = choice; break
    if trial % 50 == 0: print("trial", trial, "uncovered", left, flush=True)
print("covering found:" if best else "no covering found", best)
if best:
    import json; json.dump({"L": L, "plus": plus, "minus": {str(p): v for p, v in best.items()}}, open(f"cover69_{L}.json", "w"))
