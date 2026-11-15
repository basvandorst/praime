"""Independent recheck (no numpy sieve, no shared code): primes from gmpy2.next_prime,
features from plain integer arithmetic, Poisson-binomial by a pure-Python DP,
Poisson tail from math.  Also Lucas-Lehmer over all primes < 10^4."""
import math
import gmpy2
import sympy
from scipy.stats import kstest

EXPS = [int(x) for x in """2 3 5 7 13 17 19 31 61 89 107 127 521 607 1279 2203 2281 3217 4253 4423
9689 9941 11213 19937 21701 23209 44497 86243 110503 132049 216091 756839 859433 1257787 1398269
2976221 3021377 6972593 13466917 20996011 24036583 25964951 30402457 32582657 37156667 42643801
43112609 57885161 74207281 77232917 82589933""".split()]
SE = set(EXPS)
B = 83880077
G = math.exp(0.5772156649015329)
assert len(EXPS) == 51 and all(sympy.isprime(p) for p in EXPS)   # proved (< 2^64)


def lucas_lehmer(p):
    if p == 2:
        return True
    M = (1 << p) - 1
    s = gmpy2.mpz(4)
    for _ in range(p - 2):
        s = (s * s - 2) % M
    return s == 0


ll = [p for p in sympy.primerange(2, 10**4) if lucas_lehmer(p)]
print("Lucas-Lehmer, all primes < 10^4:", ll)
assert ll == [p for p in EXPS if p < 10**4]

# Euler: p = 3 mod 4, p > 3, q = 2p+1 prime  =>  q | 2^p - 1   (checked for p < 10^5)
for p in sympy.primerange(5, 10**5):
    if p % 4 == 3 and sympy.isprime(2 * p + 1):
        assert pow(2, p, 2 * p + 1) == 1
print("Euler divisibility checked for all p < 1e5")

NMC = set()
for k in range(30):
    for v in (2**k - 1, 2**k + 1, 4**k - 3, 4**k + 3):
        if 2 <= v < B:
            NMC.add(v)


def Pnull(p, mod4=True):
    a = 6 if (mod4 and p % 4 == 1) else 2
    return min(1.0, G * math.log2(a * p) / p)


classes = {
    "all": lambda p: True,
    "p=1 mod 4": lambda p: p % 4 == 1,
    "p=3 mod 4": lambda p: p % 4 == 3,
    "p=1 mod 8": lambda p: p % 8 == 1,
    "p=1 mod 7": lambda p: p % 7 == 1,
    "p=5 mod 12": lambda p: p % 12 == 5,
    "p=11 mod 12": lambda p: p % 12 == 11,
    "NMC form": lambda p: p in NMC,
    "near 2^k (0.1)": lambda p: min(math.log2(p) % 1, 1 - math.log2(p) % 1) < 0.1,
    "control 2p+1 prime, 3 mod 4, p>3": lambda p: p % 4 == 3 and p > 3 and gmpy2.is_prime(2 * p + 1),
    "p+2 prime": lambda p: gmpy2.is_prime(p + 2),
    "p-2 prime": lambda p: p > 2 and gmpy2.is_prime(p - 2),
    "[1e7,B)": lambda p: p >= 10**7,
    "[1e5,B)": lambda p: p >= 10**5,
}
names = list(classes)
E = {n: 0.0 for n in names}; E2 = {n: 0.0 for n in names}
Eh = {n: 0.0 for n in names}
smallP = {n: [] for n in names}       # P values for p < 1e5 in class (for exact DP)
O = {n: 0 for n in names}; Oh = {n: 0 for n in names}
# refined-null ingredients: per octave, mass of p=3 mod 4 and of the Euler-killed ones
oct_tot = [0.0] * 40; oct_kill = [0.0] * 40
# Lambda(x) for gap test
lam = 0.0; lam_at = {}
records = []  # (p, P, kill, member tuple) for p < 1e5
big3 = {n: [0.0] * 40 for n in names}    # p >= 1e5: non-killed p=3 mod 4 mass per octave, per class
bigo = {n: 0.0 for n in names}           # p >= 1e5: p != 3 mod 4 mass per class
p = 2
nprimes = 0
while p < B:
    P = Pnull(p); P2 = Pnull(p, False)
    nprimes += 1
    lam += P
    if p in SE:
        lam_at[p] = lam - P / 2
    mem = tuple(classes[n](p) for n in names)
    kill = p % 4 == 3 and p > 3 and gmpy2.is_prime(2 * p + 1)
    o = p.bit_length() - 1
    if p % 4 == 3:
        oct_tot[o] += P
        if kill:
            oct_kill[o] += P
    for n, m in zip(names, mem):
        if m:
            E[n] += P; E2[n] += P2
            if p >= 10**5:
                Eh[n] += P
            else:
                smallP[n].append(P)
            if p in SE:
                O[n] += 1
                if p >= 10**5:
                    Oh[n] += 1
    if p < 10**5:
        records.append((p, P, kill, mem))
    else:
        for n, m in zip(names, mem):
            if m and not kill:
                if p % 4 == 3:
                    big3[n][o] += P
                else:
                    bigo[n] += P
    p = int(gmpy2.next_prime(p))
print("primes < B:", nprimes)


def pb_pmf(ps, lam_rest, L=140):
    d = [1.0] + [0.0] * (L - 1)
    for q in ps:
        nd = [d[0] * (1 - q)]
        for k in range(1, L):
            nd.append(d[k] * (1 - q) + d[k - 1] * q)
        d = nd
    pois = [math.exp(-lam_rest + k * math.log(lam_rest) - math.lgamma(k + 1)) if lam_rest > 0 else float(k == 0) for k in range(L)]
    return [sum(d[i] * pois[k - i] for i in range(k + 1)) for k in range(L)]


def two_sided(pmf, o):
    lo = sum(pmf[:o + 1]); hi = sum(pmf[o:])
    return min(1.0, 2 * min(lo, hi))


def pois_two_sided(lam_, o):
    pmf = [math.exp(-lam_ + k * math.log(lam_) - math.lgamma(k + 1)) for k in range(200)]
    return two_sided(pmf, o)


print(f"\n{'class':>36} {'O':>3} {'E(a=2/6)':>9} {'E(a=2)':>8} {'p exact':>8} | {'Oh':>3} {'Eh':>6} {'ph':>6}")
for n in names:
    pe = two_sided(pb_pmf(smallP[n], Eh[n]), O[n])
    ph = pois_two_sided(Eh[n], Oh[n]) if Eh[n] > 0 else 1.0
    print(f"{n:>36} {O[n]:>3} {E[n]:>9.2f} {E2[n]:>8.2f} {pe:>8.4f} | {Oh[n]:>3} {Eh[n]:>6.2f} {ph:>6.3f}")

# refined null (Euler kills zeroed, p=3 mod 4 rescaled per octave)
f = [(oct_tot[i] / (oct_tot[i] - oct_kill[i])) if oct_tot[i] > oct_kill[i] else 1.0 for i in range(40)]
print("\nrefined null: removed mass", round(sum(oct_kill), 3))
for n in ("p=11 mod 12", "p+2 prime", "p-2 prime", "p=3 mod 4", "all"):
    j = names.index(n)
    Er_small = [0.0 if k else (min(1.0, P * f[q.bit_length() - 1]) if q % 4 == 3 else P)
                for (q, P, k, mem) in records if mem[j]]
    Er_h = bigo[n] + sum(f[i] * big3[n][i] for i in range(40))
    Er = sum(Er_small) + Er_h
    print(f"  {n}: O {O[n]}, E_ref {Er:.2f}, p_ref {two_sided(pb_pmf(Er_small, Er_h), O[n]):.3f}")

# gap test on exponents >= 1000 and >= 1e5 (Lambda-gaps vs Exp(1))
for thr in (1000, 10**5):
    u = [lam_at[p] for p in EXPS if p >= thr]
    g = [b - a for a, b in zip(u, u[1:])]
    kt = kstest(g, "expon")
    m = sum(g) / len(g); sd = math.sqrt(sum((x - m) ** 2 for x in g) / len(g))
    print(f"gaps, exponents >= {thr}: n={len(g)} D={kt.statistic:.3f} p={kt.pvalue:.3f} mean {m:.3f} CV {sd/m:.3f}")
