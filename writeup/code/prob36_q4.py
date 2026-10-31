# Problem 36 Q4: k with 3 | k such that k*2^n + 1 and k*2^n - 1 are both composite for every even n (and for every odd n).
# Even n = 2m: N = k*4^m +- 1. Odd n = 2m+1: N = (2k)*4^m +- 1. Covering systems in base 4, then CRT.
import itertools, random
from sympy import isprime, n_order, divisors
from sympy.ntheory.modular import crt
def plus_cover():
    P = {5: 2, 7: 3, 17: 4, 13: 6, 241: 12}
    ps = list(P)
    for res in itertools.product(*[range(P[p]) for p in ps]):
        if all(any(m % P[p] == r for p, r in zip(ps, res)) for m in range(12)):
            return {p: (P[p], r) for p, r in zip(ps, res)}
def minus_cover(excl, L=720, seed=1):
    avail = []
    for mod in divisors(L):
        if mod < 2: continue
        cnt = 0
        for t in range(1, 400000):
            p = t * mod + 1
            if p > 10**8 or cnt >= 6: break
            if p in excl or p == 3 or not isprime(p): continue
            if pow(4, mod, p) == 1 and n_order(4, p) == mod: avail.append((mod, p)); cnt += 1
    random.seed(seed)
    for trial in range(400):
        covered = bytearray(L); left = L; ch = {}; used = set()
        random.shuffle(avail)
        while left:
            best = None; bg = 0
            for mod, p in avail:
                if p in used: continue
                for r in range(mod):
                    g = sum(1 for n in range(r, L, mod) if not covered[n])
                    if g > bg: bg, best = g, (mod, p, r)
            if not best: break
            mod, p, r = best; used.add(p); ch[p] = (mod, r)
            for n in range(r, L, mod):
                if not covered[n]: covered[n] = 1; left -= 1
        if left == 0: return ch, L
    return None, L
plus = plus_cover(); minus, L = minus_cover(set(plus))
print("+ cover (base 4):", plus); print("- cover (base 4, mod %d): %d primes" % (L, len(minus)))
def solve(mult):            # N = (mult*k) * 4^m +- 1 covered; mult = 1 (even n) or 2 (odd n)
    mods, res = [3, 2], [0, 1]                       # 3 | k, k odd
    for p, (o, a) in plus.items():                   # p | mult*k*4^a + 1  =>  k = -(mult*4^a)^-1
        mods.append(p); res.append((-pow(mult * pow(4, a, p), -1, p)) % p)
    for q, (o, b) in minus.items():                  # q | mult*k*4^b - 1  =>  k = (mult*4^b)^-1
        mods.append(q); res.append(pow(mult * pow(4, b, q), -1, q))
    k, M = crt(mods, res); return int(k)
for mult, name in ((1, "even"), (2, "odd")):
    k = solve(mult)
    for m in range(12):
        assert any(m % o == a and (mult * k * pow(4, m, p) + 1) % p == 0 for p, (o, a) in plus.items())
    for m in range(L):
        assert any(m % o == b and (mult * k * pow(4, m, q) - 1) % q == 0 for q, (o, b) in minus.items())
    assert k % 3 == 0 and k > max(list(plus) + list(minus))
    print(f"{name} n: k = {k}  ({len(str(k))} digits, 3 | k): k*2^n +- 1 composite for every {name} n")
