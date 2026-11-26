# Referee of the merged note: Lemma 3 (squarefree level, prime 2 included), own code.
import random, math, itertools
from sympy import primerange, isprime, factorint

def kron2(a):
    # Kronecker (a/2)
    if a % 2 == 0: return 0
    return 1 if a % 8 in (1, 7) else -1

def jacobi(a, n):
    # n odd positive
    assert n > 0 and n % 2 == 1
    a %= n; res = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5): res = -res
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3: res = -res
        a %= n
    return res if n == 1 else 0

def kron(a, n):
    # Kronecker symbol (a/n), n >= 1
    res = 1
    while n % 2 == 0:
        n //= 2; res *= kron2(a)
        if res == 0: return 0
    return res * jacobi(a, n)

# cross-check own Kronecker against pari
import cypari2
pari = cypari2.Pari()
rng = random.Random(7)
for _ in range(3000):
    a = rng.randint(-10**6, 10**6); n = rng.randint(1, 5000)
    assert kron(a, n) == int(pari.kronecker(a, n)), (a, n)
print("own kronecker = pari on 3000 random pairs")

def fund(Dneg):
    # Dneg < 0, Dneg = 0,1 mod 4: return (Delta0, f0)
    f = factorint(-Dneg)
    core = -1; sq = 1
    for l, e in f.items():
        if e % 2: core *= l
        sq *= l ** (e // 2)
    # Dneg = core * sq^2
    if core % 4 == 1:
        return core, sq
    # core = 2,3 mod 4: Delta0 = 4 core, f0 = sq/2
    assert sq % 2 == 0
    return 4 * core, sq // 2

def Mset(D):
    # m >= 1 with m^2 | D and -D/m^2 = 0,1 mod 4
    f = factorint(D)
    base = 1
    for l, e in f.items(): base *= l ** (e // 2)
    out = []
    for m in range(1, base + 1):
        if base % m == 0 and D % (m * m) == 0 and (-D // (m * m)) % 4 in (0, 1):
            out.append(m)
    return out

levels = []
for q in [7, 11, 13, 101, 1009]:
    for c in [2, 6, 30]:
        N = c * q
        if all(e == 1 for e in factorint(N).values()):
            levels.append(N)
levels += [2 * 10007, 6 * 3343, 2 * 3 * 5 * 7 * 11 * 13, 2 * 3 * 5 * 7 * 11 * 13 * 17]
stats = dict(cases=0, v2={}, bad_div=0, bad_M=0, bad_T=0, Tcases=0, classes=set())
for N in levels:
    odd = [l for l in factorint(N)]
    for r in range(0, 5):
        lo = r * r * N // 4 + 1
        ps = [p for p in primerange(max(lo, 3), lo + 1500) if N % p][:60]
        for p in ps:
            u = 4 * p - r * r * N
            if u <= 0: continue
            D = N * u
            d0, f0 = fund(-D)
            stats['cases'] += 1
            if d0 % N != 0: stats['bad_div'] += 1
            if math.gcd(f0, N) != 1: stats['bad_div'] += 1
            if N % 2 == 0:
                v = 0; t = -d0
                while t % 2 == 0: t //= 2; v += 1
                stats['v2'][(r % 2, v)] = stats['v2'].get((r % 2, v), 0) + 1
                stats['classes'].add((r % 2, p % 16))
            Ms = Mset(D)
            # M(D) must equal the divisors of f0
            assert Ms == [m for m in range(1, f0 + 1) if f0 % m == 0], (N, r, p)
            if any(math.gcd(m, N) > 1 for m in Ms): stats['bad_M'] += 1
            # m = 2 explicitly for even N
            if N % 2 == 0:
                assert not ((-D // 4) % 4 in (0, 1) and D % 4 == 0), (N, r, p, (-D // 4) % 4)
            # T(m, n; p) = 0 when gcd(n, N) > 1
            for m in Ms:
                Dm = -D // (m * m)
                for _ in range(20):
                    qn = rng.choice(odd)
                    n = qn * rng.randint(1, 400)
                    stats['Tcases'] += 1
                    if kron(Dm, n) != 0: stats['bad_T'] += 1
print("Lemma 3 cases", stats['cases'], "N | Delta0 and f0 prime to N violated:", stats['bad_div'],
      " M(D) not prime to N:", stats['bad_M'], " T != 0 with (n,N)>1:", stats['bad_T'], "of", stats['Tcases'])
print("even N: (r mod 2, v2(Delta0)) counts:", stats['v2'])
print("even N: (r mod 2, p mod 16) classes covered:", len(stats['classes']), "of 16")
# -D/4 mod 4 for even N, by residues: N = 2N', r odd/even, p odd
seen = set()
for Np in range(1, 64, 2):
    for r in range(8):
        for p in range(1, 64, 2):
            N = 2 * Np; D = N * (4 * p - r * r * N)
            seen.add((r % 2, (-D // 4) % 4, D % 16 in (4, 8, 12)))
print("even N, all odd N' and p mod 64, r mod 8: (r mod 2, -D/4 mod 4, v2(D) in {2,3}):", sorted(seen))
