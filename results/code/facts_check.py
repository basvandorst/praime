# Checks for FACTS.md that run in about a minute (the one-pass statistics to 1.4e11 come from facts.c).
# Run: .venv/bin/python results/code/facts_check.py
# Independent of facts.c: its own numpy sieve to 2e8, sympy (exact below 2^64), gmpy2 (BPSW above).
import math
from fractions import Fraction
import numpy as np
import gmpy2
from sympy import isprime, nextprime, primerange, factorint, legendre_symbol, fibonacci, primitive_root, n_order, mobius

def sieve(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, math.isqrt(n) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]

def consecutive(ps):
    return all(not any(gmpy2.is_prime(m, 30) for m in range(a + 2, b, 2)) for a, b in zip(ps, ps[1:]))

# --- section 2 and 3: counts, gaps, runs below 2e8 -----------------------------------------------
P = sieve(200_000_000)
pi = {10**k: int(np.searchsorted(P, 10**k, side='right')) for k in range(1, 9)}
assert list(pi.values()) == [4, 25, 168, 1229, 9592, 78498, 664579, 5761455]
g = np.diff(P)
tw = {x: int(np.sum(g[:np.searchsorted(P, x, side='right') - 1] == 2)) for x in pi}
assert list(tw.values()) == [2, 8, 35, 205, 1224, 8169, 58980, 440312]
rec, m = [], 0                                      # maximal gaps
for i in np.nonzero(g > np.maximum.accumulate(np.concatenate(([0], g[:-1]))))[0]:
    if g[i] > m: m = int(g[i]); rec.append((m, int(P[i])))
assert rec[:12] == [(1, 2), (2, 3), (4, 7), (6, 23), (8, 89), (14, 113), (18, 523), (20, 887), (22, 1129), (34, 1327), (36, 9551), (44, 15683)]
assert rec[-5:] == [(210, 20831323), (220, 47326693), (222, 122164747), (234, 189695659), (248, 191912783)]
hist = np.bincount(g[1:])                          # gaps between odd primes below 2e8
assert int(np.argmax(hist)) == 6                   # 6 is the most common gap
r6, q6 = P[2:-1] % 6, g[2:] % 6                     # p > 3: gap = 2 (mod 6) iff p = 5, gap = 4 (mod 6) iff p = 1 (mod 6)
assert np.all((q6 == 0) | ((q6 == 2) & (r6 == 5)) | ((q6 == 4) & (r6 == 1)))
def first_runs(diffs_order):                       # first start of a run of R consecutive primes on a degree-d polynomial
    d = diffs_order; D = np.diff(P[1:].astype(np.int64), n=d + 1); z = (D == 0).astype(np.int8)
    out = {}
    zi = np.concatenate(([0], z, [0])); starts = np.nonzero(np.diff(zi) == 1)[0]; ends = np.nonzero(np.diff(zi) == -1)[0]
    for s, e in zip(starts, ends):
        R = (e - s) + d + 1
        for r in range(d + 2, R + 1):
            out.setdefault(r, int(P[1 + s]))
    return out
cp = first_runs(1)                                 # CPAP: equal gaps
assert cp[4] == 251 and cp[5] == 9843019 and cp[6] == 121174811 and max(cp) == 6
qd = first_runs(2)
assert qd[5] == 347 and qd[6] == 2903 and qd[7] == 15373 and qd[8] == 128981
# Euler-Rodriguez chains (exactly k primes A, A+2, A+6, ...), least A for k = 5..9 (puzzle 151)
er = {}
for i in np.nonzero(g == 2)[0]:                    # a chain starts with a twin pair
    L = 1
    while i + L - 1 < len(g) and g[i + L - 1] == 2 * L: L += 1
    er.setdefault(L, int(P[i]))
assert [er[k] for k in range(5, 10)] == [347, 13901, 665111, 128981, 113575727]
# Dusart (2010): a prime in (x, x(1 + 1/(2 ln^2 x))] for x >= 3275; sharp at 3274 (next prime after 3271 is 3299)
f = lambda x: x * (1 + 1 / (2 * math.log(x) ** 2))
assert f(3274) < 3299 <= f(3275)
i0 = int(np.searchsorted(P, 3299))
Pf = P[i0:-1].astype(np.float64)
assert np.all(g[i0:] <= Pf / (2 * np.log(Pf) ** 2))
MAXGAPS = [(250, 387096133), (282, 436273009), (288, 1294268491), (292, 1453168141), (320, 2300942549), (336, 3842610773),
           (354, 4302407359), (382, 10726904659), (384, 20678048297), (394, 22367084959), (456, 25056082087),
           (464, 42652618343), (468, 127976334671)]                       # maximal gaps 2e8 .. 1.4e11, from facts.c
assert all(gap <= q / (2 * math.log(q) ** 2) and isprime(q) and isprime(q + gap) and nextprime(q) == q + gap for gap, q in MAXGAPS)
# Schoenfeld (1976): a prime in (x, x(1 + 1/16597)] for x > 2010759.9; forced by the gap 148 after 2010733
assert nextprime(2010733) == 2010881 and 2010759 * (1 + 1 / 16597) < 2010881 <= 2010760 * (1 + 1 / 16597)
k0 = int(np.searchsorted(P, 2010881))
assert np.all(g[k0:] <= P[k0:-1].astype(np.float64) / 16597) and all(gap <= q / 16597 for gap, q in MAXGAPS)
# Legendre/Andrica range: g < 1.32 sqrt(p) for p >= 11 (largest ratio at 113)
j0 = int(np.searchsorted(P, 11))
ratio = g[j0:] / np.sqrt(P[j0:-1].astype(np.float64))
assert ratio.max() < 1.32 and int(P[j0 + int(np.argmax(ratio))]) == 113
assert all(gap < 1.32 * math.sqrt(q) for gap, q in MAXGAPS)
# Capelle (conjecture 53): floor(pq/r) = p + q - r for consecutive r < q < p, except for 15 values of r
r_, q_, p_ = P[:-2].astype(np.int64), P[1:-1].astype(np.int64), P[2:].astype(np.int64)
assert r_[(p_ - r_) * (q_ - r_) >= r_].tolist() == [2, 3, 5, 7, 11, 13, 19, 23, 31, 47, 83, 89, 113, 199, 1327]
# Chebyshev's bias mod 4: the first prime where pi(x;4,1) > pi(x;4,3)
d4 = np.cumsum(np.where(P[1:] % 4 == 3, 1, -1))
assert int(P[1:][np.argmax(d4 < 0)]) == 26861
# Lemke Oliver-Soundararajan: after a prime ending in 1, the next prime ends in 1 least often (first 10^7 primes)
a, b = P[3:10**7 - 1] % 10, P[4:10**7] % 10
row1 = [int(np.sum((a == 1) & (b == c))) for c in (1, 3, 7, 9)]
assert row1[0] == min(row1)

# --- admissible k-tuples: least diameter s(k) ----------------------------------------------------
def least_diameter(k, start):
    ps = list(primerange(2, k + 1))
    def add(used, h):
        new = []
        for j, p in enumerate(ps):
            mm = used[j] | (1 << (h % p))
            if mm == (1 << p) - 1: return None
            new.append(mm)
        return new
    D = start
    while True:
        u = add(add([0] * len(ps), 0), D)
        def dfs(h, c, u):
            if c == 0: return True
            for x in range(h + 2, D - 2 * (c - 1), 2):
                v = add(u, x)
                if v is not None and dfs(x, c - 1, v): return True
            return False
        if u is not None and dfs(0, k - 2, u): return D
        D += 2
s, prev = [], 2
for k in range(2, 20):
    prev = least_diameter(k, prev); s.append(prev)
assert s == [2, 6, 8, 12, 16, 20, 26, 30, 32, 36, 42, 48, 50, 56, 60, 66, 70, 76]
assert int(P[np.searchsorted(P, 113)] - 37) == 76 and len(P[(P >= 37) & (P <= 113)]) == 19    # primes 37..113: least width
assert len(P[(P >= 17) & (P <= 73)]) == 15 and 73 - 17 == 56

# --- section 4: record chains and progressions -------------------------------------------------
for k, A in [(10, 2426256797), (11, 137168442221), (12, 4656625081181), (13, 101951758179851), (14, 484511389338941)]:
    v = [A + x * x + x for x in range(k)]
    assert all(isprime(t) for t in v) and consecutive(v) and nextprime(v[-1]) != A + k * (k + 1)
p0 = 100996972469714247637786655587969840329509324689190041803603417758904341703348882159067229719
c10 = [p0 + 210 * j for j in range(10)]
assert all(gmpy2.is_prime(t, 30) for t in c10) and consecutive(c10)
assert gmpy2.next_prime(c10[-1]) != c10[-1] + 210 and gmpy2.prev_prime(p0) != p0 - 210
P23 = 223092870
ap = [224584605939537911 + 81292139 * P23 * n for n in range(27)]
assert all(isprime(t) for t in ap) and max(ap) < 2**64

# --- section 1: congruences, checked on the primes below 10^5 --------------------------------------
ps = [int(p) for p in P[P < 100000]]
two_sq = lambda n, k: any(math.isqrt((n - x * x) // k) ** 2 * k == n - x * x for x in range(math.isqrt(n) + 1) if n - x * x > 0 and (n - x * x) % k == 0)
assert all(p * p % 24 == 1 for p in ps if p > 3) and all(pow(p, 4, 240) == 1 for p in ps if p > 5)
assert math.gcd(*[p**2 - 1 for p in ps if p > 3]) == 24 and math.gcd(*[p**4 - 1 for p in ps if p > 5]) == 240
assert math.gcd(*[p**6 - 1 for p in ps if p > 7]) == 504 and math.gcd(*[p**8 - 1 for p in ps if p > 5]) == 480
small = [p for p in ps if p < 10000]
assert all((p % 4 == 1) == two_sq(p, 1) for p in small if p > 2)
assert all((p == 2 or p % 8 in (1, 3)) == two_sq(p, 2) for p in small)
assert all((p == 3 or p % 3 == 1) == two_sq(p, 3) for p in small)
assert all((p % 20 in (1, 9)) == two_sq(p, 5) for p in small if p not in (2, 5))
assert all((p % 7 in (1, 2, 4)) == two_sq(p, 7) for p in small if p not in (2, 7))
cubic2 = lambda p: pow(2, (p - 1) // 3, p) == 1
quart2 = lambda p: pow(2, (p - 1) // 4, p) == 1
assert all((p % 3 == 1 and cubic2(p)) == two_sq(p, 27) for p in small if p > 3)
assert all((p % 8 == 1 and quart2(p)) == two_sq(p, 64) for p in small if p > 2)
for p in small:
    if p % 4 == 1:                                       # Gauss: C((p-1)/2, (p-1)/4) = 2a (mod p), p = a^2 + b^2, a = 1 (mod 4)
        a_ = next(x for x in range(1, math.isqrt(p) + 1, 2) if math.isqrt(p - x * x) ** 2 == p - x * x)
        a_ = a_ if a_ % 4 == 1 else -a_
        assert math.comb((p - 1) // 2, (p - 1) // 4) % p == (2 * a_) % p
        assert pow(math.factorial((p - 1) // 2), 2, p) == p - 1          # ((p-1)/2)!^2 = -1 (mod p)
assert all(pow(math.factorial((p - 1) // 2), 2, p) == 1 for p in small if p % 4 == 3)
assert all(math.comb(2 * p - 1, p - 1) % p**3 == 1 for p in ps if 5 <= p < 2000)
assert all(fibonacci(p - legendre_symbol(p % 5, 5)) % p == 0 and fibonacci(p) % p == legendre_symbol(p % 5, 5) % p for p in small if p not in (2, 5))
assert all((legendre_symbol(2, p) == 1) == (p % 8 in (1, 7)) for p in ps if p > 2)
assert all((legendre_symbol(3, p) == 1) == (p % 12 in (1, 11)) for p in ps if p > 3)
assert all((legendre_symbol(5, p) == 1) == (p % 5 in (1, 4)) for p in ps if p > 5)
assert all((legendre_symbol(p - 3, p) == 1) == (p % 3 == 1) for p in ps if p > 3)
assert all(((4 * (math.factorial(n - 1) + 1) + n) % (n * (n + 2)) == 0) == (isprime(n) and isprime(n + 2)) for n in range(3, 1500, 2))   # Clement
assert [p for p in primerange(3, 10**7) if pow(2, p - 1, p * p) == 1] == [1093, 3511]                         # Wieferich
wil, f = [], 1
for n in range(2, 1001):
    f *= n - 1                                                                                                 # f = (n-1)!
    if isprime(n) and (f + 1) % (n * n) == 0: wil.append(n)
assert wil == [5, 13, 563]                                                                                   # Wilson primes
assert all(all(q % (2 * p) == 1 and q % 8 in (1, 7) for q in factorint(2**p - 1)) for p in primerange(3, 60))
assert all(pow(2, p, 2 * p + 1) == 1 for p in ps if p % 4 == 3 and p > 3 and isprime(2 * p + 1))
def lucas_lehmer(p):
    M, s = (1 << p) - 1, 4
    for _ in range(p - 2): s = (s * s - 2) % M
    return s == 0
assert [p for p in primerange(3, 1300) if lucas_lehmer(p)] == [3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279]
pepin = lambda n: gmpy2.powmod(3, 2**(2**n - 1), 2**(2**n) + 1) == 2**(2**n)      # 3^((F-1)/2) = -1 (mod F)
assert [n for n in range(1, 15) if pepin(n)] == [1, 2, 3, 4]
assert all(q % 2**(n + 2) == 1 for n in range(2, 7) for q in factorint(2**(2**n) + 1))   # Lucas: factors of F_n
art = [p for p in ps[1:] if n_order(2, p) == p - 1]
assert abs(len(art) / (len(ps) - 1) - 0.3739558) < 0.01                             # Artin's constant
r, prev, diffs = 7, 7, set()                                                       # Rowland: a(n) = a(n-1) + gcd(n, a(n-1))
for n in range(2, 200000):
    r = prev + math.gcd(n, prev); d = r - prev
    assert d == 1 or isprime(d); diffs.add(d); prev = r
# --- section 10: facts from the archive that we can recompute -----------------------------------
assert [p for p in primerange(5, 10**6) if isprime((p * p - 1) // 24)] == [7, 11, 13]           # (p^2-1)/24 prime only 2, 5, 7
form = {6 * a * b + sa * a + sb * b for a in range(1, 2000) for b in range(1, 2000) for sa in (1, -1) for sb in (1, -1) if 6 * a * b - a - b <= 10**4}
assert all((isprime(6 * m - 1) and isprime(6 * m + 1)) == (m not in form) for m in range(1, 10**4))  # Perry / Balestrieri
cover = [3, 5, 7, 13, 19, 37, 73]                                                                 # Selfridge: 78557*2^n + 1
assert all(any((78557 * pow(2, n, q) + 1) % q == 0 for q in cover) for n in range(1, 145))
for p in primerange(5, 2000):                                                                      # primitive roots: product 1, sum mu(p-1)
    roots = [g for g in range(1, p) if n_order(g, p) == p - 1]
    assert math.prod(roots) % p == 1 and sum(roots) % p == mobius(p - 1) % p
for p in primerange(3, 400):                                                                       # Lerch: sum of Fermat quotients = Wilson quotient
    assert sum((pow(a, p - 1, p * p) - 1) // p for a in range(1, p)) % p == ((math.factorial(p - 1) + 1) // p) % p
assert [p for p in [3, 5, 11, 29] if all(isprime(2 * n * n + p) for n in range(p))] == [3, 5, 11, 29]
assert [p for p in primerange(3, 200) if all(isprime(2 * n * n + p) for n in range(p))] == [3, 5, 11, 29]   # Euler 2n^2 + p
def ram(n, X=2000):                                                                                # Ramanujan primes
    pc = np.cumsum(np.isin(np.arange(X + 1), P[P <= X]))
    return max(x for x in range(1, X + 1) if pc[x] - pc[x // 2] < n) + 1
assert [ram(n) for n in range(1, 6)] == [2, 11, 17, 29, 41]
assert all((int(str(p)[::-1]) - p) % 18 == 0 for p in ps if isprime(int(str(p)[::-1])))          # emirp pairs differ by 18k
assert all(fibonacci(p - 1) % p == 0 for p in small if p % 5 in (1, 4)) and all(fibonacci(p + 1) % p == 0 for p in small if p % 5 in (2, 3))
assert [p for p in ps if p > 5 and n_order(10, p) == p - 1][:6] == [7, 17, 19, 23, 29, 47]          # full reptend primes
assert isprime(18691113008663) and isprime(18691113008663 + 666) and nextprime(18691113008663) == 18691113008663 + 666   # first gap 666
is2 = (g == 2)                                                                                      # consecutive twin pairs (puzzle 122)
def first_twin_run(k):                                                                             # k twin pairs, no other primes between:
    ok = is2[:len(is2) - 2 * (k - 1)].copy()                                                       # gaps g[s], g[s+2], ..., g[s+2(k-1)] all 2
    for j in range(1, k): ok &= is2[2 * j:len(is2) - 2 * (k - 1) + 2 * j]
    return int(P[np.argmax(ok)])
assert [first_twin_run(k) for k in (3, 4, 5)] == [5, 9419, 909287]
d20 = next(i for i in range(len(g) - 20) if len(set(g[i:i + 20].tolist())) == 20)                  # 20 distinct gaps in a row
assert int(P[d20]) == 60155567
assert [int(P[n - 1]) for n in range(1, 200000) if n * math.log(int(P[n]) / int(P[n - 1])) > math.log(n)] == [2, 3, 7, 113, 1327]
eu, pr = [], 1
for q in primerange(2, 20):
    pr *= q; eu.append(isprime(pr + 1))
assert eu[:5] == [True] * 5 and not eu[5] and 2 * 3 * 5 * 7 * 11 * 13 + 1 == 59 * 509
odd_comp = [n for n in range(9, 10**6, 2) if not isprime(n)]                                       # Goldbach's other conjecture
assert [n for n in odd_comp if not any(isprime(n - 2 * k * k) for k in range(1, math.isqrt(n // 2) + 1))] == [5777, 5993]
row = [int(x) for x in P[:2000]]                                                                    # Gilbreath on the first 2000 primes
for _ in range(1999):
    row = [abs(b - a) for a, b in zip(row, row[1:])]
    assert row[0] == 1 or len(row) < 2
    if len(row) < 2: break
# --- section 5: formulas -----------------------------------------------------------------------
good = [A for A in range(2, 10001) if all(isprime(x * x + x + A) for x in range(A - 1))]
assert good == [2, 3, 5, 11, 17, 41]                                               # Rabinowitsch / Heegner
assert all(legendre_symbol(-163 % p, p) == -1 for p in primerange(3, 41))         # no prime < 41 divides n^2 + n + 41
ms = [2]
for _ in range(5): ms.append(int(gmpy2.next_prime(ms[-1] ** 3)))
assert ms[:5] == [2, 11, 1361, 2521008887, 16022236204009818131831320183]
gmpy2.get_context().precision = 2000
assert str(gmpy2.root(gmpy2.mpfr(ms[5]), 729))[:22] == "1.30637788386308069046"
fr, prod = Fraction(0), 1
for p in primerange(2, 600): fr += Fraction(p - 1, prod); prod *= p
assert str(fr.numerator * 10**30 // fr.denominator)[:25] == "2920050977316134712092562"
x, got = fr, []
for _ in range(60):
    q = x.numerator // x.denominator; got.append(q); x = q * (x - q + 1)
assert got == list(primerange(2, 300))[:60]
def willans(n):
    s, tot, fct = 0, 0, 1
    for i in range(1, 2**n + 1):
        if i > 1: fct *= i - 1
        tot += 1 if (fct + 1) % i == 0 else 0
        s += 1 if tot <= n else 0                       # floor((n/tot)^(1/n)) is 1 iff tot <= n, else 0
    return 1 + s
assert [willans(n) for n in range(1, 9)] == [2, 3, 5, 7, 11, 13, 17, 19]
# --- constants -----------------------------------------------------------------------------------
Q = P[P < 10**7].astype(np.float64)
C2 = float(np.prod(1 - 1 / (Q[1:] - 1) ** 2)); assert abs(C2 - 0.6601618158) < 1e-7
ART = float(np.prod(1 - 1 / (Q * (Q - 1)))); assert abs(ART - 0.3739558136) < 1e-7
MM = float(np.sum(1 / P.astype(np.float64)) - math.log(math.log(2e8))); assert abs(MM - 0.2614972128) < 1e-4
print("all FACTS.md checks passed")
