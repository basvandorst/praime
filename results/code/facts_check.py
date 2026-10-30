# Checks for FACTS.md, about 3 minutes (the one-pass statistics to 1.4e11 come from the C programs in this folder).
# Run: .venv/bin/python results/code/facts_check.py
# Independent of the C programs: its own numpy sieve to 2e8, sympy (exact below 2^64), gmpy2 (BPSW above).
# Every check starts with a comment "# NNNN: what it checks", where NNNN is the fact number in FACTS.md.
import math
from fractions import Fraction
from itertools import permutations
import numpy as np
import gmpy2
import mpmath
from sympy import isprime, nextprime, primerange, factorint, legendre_symbol, fibonacci, lucas, primitive_root, n_order, mobius
from sympy import cyclotomic_poly, primefactors, totient

def sieve(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, math.isqrt(n) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]

def consecutive(ps):
    return all(not any(gmpy2.is_prime(m, 30) for m in range(a + 2, b, 2)) for a, b in zip(ps, ps[1:]))

P = sieve(200_000_000)
g = np.diff(P)
isp = np.zeros(200_000_001, dtype=bool); isp[P] = True

def at(a): return np.where(a <= 200_000_000, isp[np.minimum(a, 200_000_000)], False)

# --- sections 2 and 3: counts, gaps, runs below 2e8 ---------------------------------------------------

# 2011: pi(10^k) for k = 1..8
pi = {10**k: int(np.searchsorted(P, 10**k, side='right')) for k in range(1, 9)}
assert list(pi.values()) == [4, 25, 168, 1229, 9592, 78498, 664579, 5761455]

# 3014: twin pairs below 10^k
tw = {x: int(np.sum(g[:np.searchsorted(P, x, side='right') - 1] == 2)) for x in pi}
assert list(tw.values()) == [2, 8, 35, 205, 1224, 8169, 58980, 440312]

# 3011: maximal gaps below 2e8
rec, m = [], 0
for i in np.nonzero(g > np.maximum.accumulate(np.concatenate(([0], g[:-1]))))[0]:
    if g[i] > m: m = int(g[i]); rec.append((m, int(P[i])))
assert rec[:12] == [(1, 2), (2, 3), (4, 7), (6, 23), (8, 89), (14, 113), (18, 523), (20, 887), (22, 1129), (34, 1327), (36, 9551), (44, 15683)]
assert rec[-5:] == [(210, 20831323), (220, 47326693), (222, 122164747), (234, 189695659), (248, 191912783)]

# 3013: 6 is the most common gap between odd primes below 2e8
hist = np.bincount(g[1:])
assert int(np.argmax(hist)) == 6

# 1108: for p > 3, gap = 2 (mod 6) iff p = 5 (mod 6), gap = 4 (mod 6) iff p = 1 (mod 6)
r6, q6 = P[2:-1] % 6, g[2:] % 6
assert np.all((q6 == 0) | ((q6 == 2) & (r6 == 5)) | ((q6 == 4) & (r6 == 1)))

# 4007, 4011: first runs of consecutive primes in arithmetic progression and on a quadratic
def first_runs(diffs_order):                       # first start of a run of R consecutive primes on a degree-d polynomial
    d = diffs_order; D = np.diff(P[1:].astype(np.int64), n=d + 1); z = (D == 0).astype(np.int8)
    out = {}
    zi = np.concatenate(([0], z, [0])); starts = np.nonzero(np.diff(zi) == 1)[0]; ends = np.nonzero(np.diff(zi) == -1)[0]
    for s, e in zip(starts, ends):
        R = (e - s) + d + 1
        for r in range(d + 2, R + 1):
            out.setdefault(r, int(P[1 + s]))
    return out
cp = first_runs(1)
assert cp[4] == 251 and cp[5] == 9843019 and cp[6] == 121174811 and max(cp) == 6
qd = first_runs(2)
assert qd[5] == 347 and qd[6] == 2903 and qd[7] == 15373 and qd[8] == 128981

# 4009: Euler-Rodriguez chains of exactly k primes A, A+2, A+6, ..., least A for k = 5..9
er = {}
for i in np.nonzero(g == 2)[0]:                    # a chain starts with a twin pair
    L = 1
    while i + L - 1 < len(g) and g[i + L - 1] == 2 * L: L += 1
    er.setdefault(L, int(P[i]))
assert [er[k] for k in range(5, 10)] == [347, 13901, 665111, 128981, 113575727]

# 3003: Dusart 2010, a prime in (x, x(1 + 1/(2 ln^2 x))] for x >= 3275, sharp at 3274
f = lambda x: x * (1 + 1 / (2 * math.log(x) ** 2))
assert f(3274) < 3299 <= f(3275)
i0 = int(np.searchsorted(P, 3299))
Pf = P[i0:-1].astype(np.float64)
assert np.all(g[i0:] <= Pf / (2 * np.log(Pf) ** 2))

# 3011, 3003: maximal gaps from 2e8 to 1.4e11 (from facts.c) are prime gaps and obey Dusart's bound
MAXGAPS = [(250, 387096133), (282, 436273009), (288, 1294268491), (292, 1453168141), (320, 2300942549), (336, 3842610773),
           (354, 4302407359), (382, 10726904659), (384, 20678048297), (394, 22367084959), (456, 25056082087),
           (464, 42652618343), (468, 127976334671)]
assert all(gap <= q / (2 * math.log(q) ** 2) and isprime(q) and isprime(q + gap) and nextprime(q) == q + gap for gap, q in MAXGAPS)

# 10201: Schoenfeld 1976, a prime in (x, x(1 + 1/16597)] for x > 2010759.9, forced by the gap 148 after 2010733
assert nextprime(2010733) == 2010881 and 2010759 * (1 + 1 / 16597) < 2010881 <= 2010760 * (1 + 1 / 16597)
k0 = int(np.searchsorted(P, 2010881))
assert np.all(g[k0:] <= P[k0:-1].astype(np.float64) / 16597) and all(gap <= q / 16597 for gap, q in MAXGAPS)

# 3006: g < 1.32 sqrt(p) for p >= 11, largest ratio at 113
j0 = int(np.searchsorted(P, 11))
ratio = g[j0:] / np.sqrt(P[j0:-1].astype(np.float64))
assert ratio.max() < 1.32 and int(P[j0 + int(np.argmax(ratio))]) == 113
assert all(gap < 1.32 * math.sqrt(q) for gap, q in MAXGAPS)

# 3017: Capelle, floor(pq/r) = p + q - r for consecutive r < q < p except for 15 values of r
r_, q_, p_ = P[:-2].astype(np.int64), P[1:-1].astype(np.int64), P[2:].astype(np.int64)
assert r_[(p_ - r_) * (q_ - r_) >= r_].tolist() == [2, 3, 5, 7, 11, 13, 19, 23, 31, 47, 83, 89, 113, 199, 1327]

# 1602: Chebyshev's bias, the first prime where pi(x;4,1) > pi(x;4,3) is 26861
d4 = np.cumsum(np.where(P[1:] % 4 == 3, 1, -1))
assert int(P[1:][np.argmax(d4 < 0)]) == 26861

# 1604: after a prime ending in 1, the next prime ends in 1 least often (first 10^7 primes)
a, b = P[3:10**7 - 1] % 10, P[4:10**7] % 10
row1 = [int(np.sum((a == 1) & (b == c))) for c in (1, 3, 7, 9)]
assert row1[0] == min(row1)

# --- section 4: admissible tuples, chains and progressions ------------------------------------------

# 4002: least diameter s(k) of an admissible k-tuple, k = 2..19
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

# 4003: the 19 primes 37..113 and the 15 primes 17..73 have the least possible width
assert int(P[np.searchsorted(P, 113)] - 37) == 76 and len(P[(P >= 37) & (P <= 113)]) == 19
assert len(P[(P >= 17) & (P <= 73)]) == 15 and 73 - 17 == 56

# 4009: Andersen's Euler-Rodriguez chains of 10..14 consecutive primes A + x(x + 1)
for k, A in [(10, 2426256797), (11, 137168442221), (12, 4656625081181), (13, 101951758179851), (14, 484511389338941)]:
    v = [A + x * x + x for x in range(k)]
    assert all(isprime(t) for t in v) and consecutive(v) and nextprime(v[-1]) != A + k * (k + 1)

# 4007: the 1998 CPAP-10 (93 digits) is ten consecutive PRPs and not longer
p0 = 100996972469714247637786655587969840329509324689190041803603417758904341703348882159067229719
c10 = [p0 + 210 * j for j in range(10)]
assert all(gmpy2.is_prime(t, 30) for t in c10) and consecutive(c10)
assert gmpy2.next_prime(c10[-1]) != c10[-1] + 210 and gmpy2.prev_prime(p0) != p0 - 210

# 4008: PrimeGrid's AP-27, all terms below 2^64 and prime
P23 = 223092870
ap = [224584605939537911 + 81292139 * P23 * n for n in range(27)]
assert all(isprime(t) for t in ap) and max(ap) < 2**64

# --- section 1: congruences, checked on the primes below 10^5 ----------------------------------------
ps = [int(p) for p in P[P < 100000]]
small = [p for p in ps if p < 10000]
two_sq = lambda n, k: any(math.isqrt((n - x * x) // k) ** 2 * k == n - x * x for x in range(math.isqrt(n) + 1) if n - x * x > 0 and (n - x * x) % k == 0)

# 1104, 1105: p^2 = 1 (mod 24), and the largest moduli for p^2, p^4, p^6, p^8
assert all(p * p % 24 == 1 for p in ps if p > 3) and all(pow(p, 4, 240) == 1 for p in ps if p > 5)
assert math.gcd(*[p**2 - 1 for p in ps if p > 3]) == 24 and math.gcd(*[p**4 - 1 for p in ps if p > 5]) == 240
assert math.gcd(*[p**6 - 1 for p in ps if p > 7]) == 504 and math.gcd(*[p**8 - 1 for p in ps if p > 5]) == 480

# 1306: p = a^2 + b^2 iff p = 1 (mod 4)
assert all((p % 4 == 1) == two_sq(p, 1) for p in small if p > 2)

# 1307: primes x^2 + n y^2 for n = 2, 3, 5, 7, 27, 64
assert all((p == 2 or p % 8 in (1, 3)) == two_sq(p, 2) for p in small)
assert all((p == 3 or p % 3 == 1) == two_sq(p, 3) for p in small)
assert all((p % 20 in (1, 9)) == two_sq(p, 5) for p in small if p not in (2, 5))
assert all((p % 7 in (1, 2, 4)) == two_sq(p, 7) for p in small if p not in (2, 7))
cubic2 = lambda p: pow(2, (p - 1) // 3, p) == 1
quart2 = lambda p: pow(2, (p - 1) // 4, p) == 1
assert all((p % 3 == 1 and cubic2(p)) == two_sq(p, 27) for p in small if p > 3)
assert all((p % 8 == 1 and quart2(p)) == two_sq(p, 64) for p in small if p > 2)

# 1306, 1203: Gauss's binomial formula for a in p = a^2 + b^2, and ((p-1)/2)!^2 = -1 (mod p) for p = 1 (mod 4)
for p in small:
    if p % 4 == 1:
        a_ = next(x for x in range(1, math.isqrt(p) + 1, 2) if math.isqrt(p - x * x) ** 2 == p - x * x)
        a_ = a_ if a_ % 4 == 1 else -a_
        assert math.comb((p - 1) // 2, (p - 1) // 4) % p == (2 * a_) % p
        assert pow(math.factorial((p - 1) // 2), 2, p) == p - 1

# 1203: ((p-1)/2)!^2 = 1 (mod p) for p = 3 (mod 4)
assert all(pow(math.factorial((p - 1) // 2), 2, p) == 1 for p in small if p % 4 == 3)

# 1204: Wolstenholme, C(2p-1, p-1) = 1 (mod p^3) for 5 <= p < 2000
assert all(math.comb(2 * p - 1, p - 1) % p**3 == 1 for p in ps if 5 <= p < 2000)

# 1206: F_p = (p/5) and p | F_(p - (p/5))
assert all(fibonacci(p - legendre_symbol(p % 5, 5)) % p == 0 and fibonacci(p) % p == legendre_symbol(p % 5, 5) % p for p in small if p not in (2, 5))

# 1303: when 2, 3, 5 and -3 are squares mod p
assert all((legendre_symbol(2, p) == 1) == (p % 8 in (1, 7)) for p in ps if p > 2)
assert all((legendre_symbol(3, p) == 1) == (p % 12 in (1, 11)) for p in ps if p > 3)
assert all((legendre_symbol(5, p) == 1) == (p % 5 in (1, 4)) for p in ps if p > 5)
assert all((legendre_symbol(p - 3, p) == 1) == (p % 3 == 1) for p in ps if p > 3)

# 1209: Clement's criterion for twin primes, odd n < 1500
assert all(((4 * (math.factorial(n - 1) + 1) + n) % (n * (n + 2)) == 0) == (isprime(n) and isprime(n + 2)) for n in range(3, 1500, 2))

# 1207: the Wieferich primes below 10^7 are 1093 and 3511
assert [p for p in primerange(3, 10**7) if pow(2, p - 1, p * p) == 1] == [1093, 3511]

# 1208: the Wilson primes below 1000 are 5, 13, 563
wil, f = [], 1
for n in range(2, 1001):
    f *= n - 1                                     # f = (n-1)!
    if isprime(n) and (f + 1) % (n * n) == 0: wil.append(n)
assert wil == [5, 13, 563]

# 1503: prime factors q of 2^p - 1 are 1 (mod 2p) and +-1 (mod 8)
assert all(all(q % (2 * p) == 1 and q % 8 in (1, 7) for q in factorint(2**p - 1)) for p in primerange(3, 60))

# 1505: 2p + 1 divides 2^p - 1 when p = 3 (mod 4) and 2p + 1 is prime
assert all(pow(2, p, 2 * p + 1) == 1 for p in ps if p % 4 == 3 and p > 3 and isprime(2 * p + 1))

# 1507, 6101: the Lucas-Lehmer test finds the Mersenne exponents up to 1279
def lucas_lehmer(p):
    M, s = (1 << p) - 1, 4
    for _ in range(p - 2): s = (s * s - 2) % M
    return s == 0
assert [p for p in primerange(3, 1300) if lucas_lehmer(p)] == [3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279]

# 1502: Pepin's test, F_1..F_4 prime and F_5..F_14 composite
pepin = lambda n: gmpy2.powmod(3, 2**(2**n - 1), 2**(2**n) + 1) == 2**(2**n)
assert [n for n in range(1, 15) if pepin(n)] == [1, 2, 3, 4]

# 1504: prime factors of F_n are 1 (mod 2^(n+2)), n = 2..6
assert all(q % 2**(n + 2) == 1 for n in range(2, 7) for q in factorint(2**(2**n) + 1))

# 1402: 2 is a primitive root for about 37.4 % of the primes below 10^5
art = [p for p in ps[1:] if n_order(2, p) == p - 1]
assert abs(len(art) / (len(ps) - 1) - 0.3739558) < 0.01

# 1210: Rowland's sequence, every difference is 1 or a prime
r, prev, diffs = 7, 7, set()
for n in range(2, 200000):
    r = prev + math.gcd(n, prev); d = r - prev
    assert d == 1 or isprime(d); diffs.add(d); prev = r

# --- section 10: facts from the archive that we can recompute -----------------------------------------

# 10101: (p^2 - 1)/24 is prime only for p = 7, 11, 13
assert [p for p in primerange(5, 10**6) if isprime((p * p - 1) // 24)] == [7, 11, 13]

# 10102: 6m +- 1 are twin primes iff m is not 6ab +- a +- b
form = {6 * a * b + sa * a + sb * b for a in range(1, 2000) for b in range(1, 2000) for sa in (1, -1) for sb in (1, -1) if 6 * a * b - a - b <= 10**4}
assert all((isprime(6 * m - 1) and isprime(6 * m + 1)) == (m not in form) for m in range(1, 10**4))

# 10113: Selfridge's covering set for 78557 * 2^n + 1
cover = [3, 5, 7, 13, 19, 37, 73]
assert all(any((78557 * pow(2, n, q) + 1) % q == 0 for q in cover) for n in range(1, 145))

# 10108: the product of the primitive roots is 1 and their sum is mu(p - 1) (mod p)
for p in primerange(5, 2000):
    roots = [g for g in range(1, p) if n_order(g, p) == p - 1]
    assert math.prod(roots) % p == 1 and sum(roots) % p == mobius(p - 1) % p

# 10105: Lerch, the sum of the Fermat quotients equals the Wilson quotient (mod p)
for p in primerange(3, 400):
    assert sum((pow(a, p - 1, p * p) - 1) // p for a in range(1, p)) % p == ((math.factorial(p - 1) + 1) // p) % p

# 10602: 2n^2 + p is prime for n = 0..p-1 exactly for p = 3, 5, 11, 29
assert [p for p in [3, 5, 11, 29] if all(isprime(2 * n * n + p) for n in range(p))] == [3, 5, 11, 29]
assert [p for p in primerange(3, 200) if all(isprime(2 * n * n + p) for n in range(p))] == [3, 5, 11, 29]

# 10201: the first Ramanujan primes
def ram(n, X=2000):
    pc = np.cumsum(np.isin(np.arange(X + 1), P[P <= X]))
    return max(x for x in range(1, X + 1) if pc[x] - pc[x // 2] < n) + 1
assert [ram(n) for n in range(1, 6)] == [2, 11, 17, 29, 41]

# 10111: emirp pairs differ by a multiple of 18
assert all((int(str(p)[::-1]) - p) % 18 == 0 for p in ps if isprime(int(str(p)[::-1])))

# 10107: p | F_(p-1) iff p = +-1 (mod 5), p | F_(p+1) iff p = +-2 (mod 5)
assert all(fibonacci(p - 1) % p == 0 for p in small if p % 5 in (1, 4)) and all(fibonacci(p + 1) % p == 0 for p in small if p % 5 in (2, 3))

# 10109: the first full reptend primes
assert [p for p in ps if p > 5 and n_order(10, p) == p - 1][:6] == [7, 17, 19, 23, 29, 47]

# 10208: the first gap of 666 follows 18691113008663
assert isprime(18691113008663) and isprime(18691113008663 + 666) and nextprime(18691113008663) == 18691113008663 + 666

# 10303: least starts of 3, 4, 5 consecutive twin pairs
is2 = (g == 2)
def first_twin_run(k):                             # gaps g[s], g[s+2], ..., g[s+2(k-1)] all 2
    ok = is2[:len(is2) - 2 * (k - 1)].copy()
    for j in range(1, k): ok &= is2[2 * j:len(is2) - 2 * (k - 1) + 2 * j]
    return int(P[np.argmax(ok)])
assert [first_twin_run(k) for k in (3, 4, 5)] == [5, 9419, 909287]

# 10304: the first 20 distinct consecutive gaps start at 60155567
d20 = next(i for i in range(len(g) - 20) if len(set(g[i:i + 20].tolist())) == 20)
assert int(P[d20]) == 60155567

# 10803: the primes with (p_(n+1)/p_n)^n > n below 2.7e6
assert [int(P[n - 1]) for n in range(1, 200000) if n * math.log(int(P[n]) / int(P[n - 1])) > math.log(n)] == [2, 3, 7, 113, 1327]

# 10704: Euclid numbers 1 + p# are prime up to 11, and 30031 = 59 * 509
eu, pr = [], 1
for q in primerange(2, 20):
    pr *= q; eu.append(isprime(pr + 1))
assert eu[:5] == [True] * 5 and not eu[5] and 2 * 3 * 5 * 7 * 11 * 13 + 1 == 59 * 509

# 10704: Goldbach's other conjecture fails only at 5777 and 5993 below 10^6
odd_comp = [n for n in range(9, 10**6, 2) if not isprime(n)]
assert [n for n in odd_comp if not any(isprime(n - 2 * k * k) for k in range(1, math.isqrt(n // 2) + 1))] == [5777, 5993]

# 10209: Gilbreath's conjecture on the first 2000 primes
row = [int(x) for x in P[:2000]]
for _ in range(1999):
    row = [abs(b - a) for a, b in zip(row, row[1:])]
    assert row[0] == 1 or len(row) < 2
    if len(row) < 2: break

# --- section 5: formulas ------------------------------------------------------------------------------

# 5003: n^2 + n + A is prime for n = 0..A-2 only for A = 2, 3, 5, 11, 17, 41
good = [A for A in range(2, 10001) if all(isprime(x * x + x + A) for x in range(A - 1))]
assert good == [2, 3, 5, 11, 17, 41]

# 5002: no prime below 41 divides n^2 + n + 41
assert all(legendre_symbol(-163 % p, p) == -1 for p in primerange(3, 41))

# 5012: Mills' constant from the least Mills primes
ms = [2]
for _ in range(5): ms.append(int(gmpy2.next_prime(ms[-1] ** 3)))
assert ms[:5] == [2, 11, 1361, 2521008887, 16022236204009818131831320183]
gmpy2.get_context().precision = 2000
assert str(gmpy2.root(gmpy2.mpfr(ms[5]), 729))[:22] == "1.30637788386308069046"

# 5012: Fridman's constant and its recursion return the first 60 primes
fr, prod = Fraction(0), 1
for p in primerange(2, 600): fr += Fraction(p - 1, prod); prod *= p
assert str(fr.numerator * 10**30 // fr.denominator)[:25] == "2920050977316134712092562"
x, got = fr, []
for _ in range(60):
    q = x.numerator // x.denominator; got.append(q); x = q * (x - q + 1)
assert got == list(primerange(2, 300))[:60]

# 5012: Willans' formula for n <= 8
def willans(n):
    s, tot, fct = 0, 0, 1
    for i in range(1, 2**n + 1):
        if i > 1: fct *= i - 1
        tot += 1 if (fct + 1) % i == 0 else 0
        s += 1 if tot <= n else 0                  # floor((n/tot)^(1/n)) is 1 iff tot <= n, else 0
    return 1 + s
assert [willans(n) for n in range(1, 9)] == [2, 3, 5, 7, 11, 13, 17, 19]

# --- section 7: constants -----------------------------------------------------------------------------

# 7001: twin prime constant, Artin's constant, Meissel-Mertens constant
Q = P[P < 10**7].astype(np.float64)
C2 = float(np.prod(1 - 1 / (Q[1:] - 1) ** 2)); assert abs(C2 - 0.6601618158) < 1e-7
ART = float(np.prod(1 - 1 / (Q * (Q - 1)))); assert abs(ART - 0.3739558136) < 1e-7
MM = float(np.sum(1 / P.astype(np.float64)) - math.log(math.log(2e8))); assert abs(MM - 0.2614972128) < 1e-4

# --- first research pass (October 2026) ---------------------------------------------------------------

# 4006: cousin, sexy, Sophie Germain, triplet and quadruplet counts below 10^k
def upto(Q, k): return Q[Q <= 10**k]
tup = {name: [int(f(upto(P, k)).sum()) for k in range(1, 9)] for name, f in {
    'cousin': lambda Q: at(Q + 4), 'sexy': lambda Q: at(Q + 6), 'sg': lambda Q: at(2 * Q + 1),
    'trip': lambda Q: (at(Q + 2) & at(Q + 6)) | (at(Q + 4) & at(Q + 6)), 'quad': lambda Q: at(Q + 2) & at(Q + 6) & at(Q + 8)}.items()}
assert tup['cousin'] == [2, 9, 41, 203, 1216, 8144, 58622, 440258] and tup['sexy'] == [2, 16, 74, 411, 2447, 16386, 117207, 879908]
assert tup['sg'] == [3, 10, 37, 190, 1171, 7746, 56032, 423140] and tup['trip'] == [2, 9, 30, 112, 507, 2837, 17220, 111156]
assert tup['quad'] == [1, 2, 5, 12, 38, 166, 899, 4768]

# 6201: least Cunningham chains of >= k primes, and their residues mod 30
def cc_least(d):                                   # chains p, 2p + d, ...: least start with >= k primes
    pre = (P - d) // 2; start = ~((pre >= 2) & at(np.maximum(pre, 0)) & (2 * pre + d == P))
    cur = P[start]; st = cur.copy(); L = np.ones(len(cur), np.int64); alive = np.ones(len(cur), bool)
    while alive.any():
        cur = 2 * cur + d; ok = alive & at(cur)
        for i in np.nonzero(alive & (cur > 200_000_000))[0]: ok[i] = gmpy2.is_prime(int(cur[i]))
        alive = ok; L += alive
    return {k: int(st[np.argmax(L >= k)]) for k in range(5, 10) if (L >= k).any()}
assert cc_least(1) == {5: 2, 6: 89, 7: 1122659, 8: 19099919, 9: 85864769}
assert cc_least(-1) == {5: 1531, 6: 16651, 7: 16651, 8: 15514861}
assert all(x % 30 == 29 for x in (89, 1122659, 19099919, 85864769)) and all(x % 30 == 1 for x in (1531, 16651, 15514861))

# 6201: Armitage's chains of 20 primes (2026) are complete
for p, d in ((75949683268775452969189499, 1), (3352894966466100254851081, -1)):
    L, x = 0, p
    while gmpy2.is_prime(x): L += 1; x = 2 * x + d
    assert L == 20 and not gmpy2.is_prime((p - d) // 2)

# 1605: the class of 1 is strictly last in the mod 8 race from 23 and in the mod 12 race from 463
for m, cl, x0 in ((8, (1, 3, 5, 7), 23), (12, (1, 5, 7, 11), 463)):
    C = np.stack([np.cumsum(P % m == c) for c in cl]); i0 = int(np.searchsorted(P, x0))
    assert np.all(C[0, i0:] < C[1:, i0:].min(0)) and not np.all(C[0, i0 - 1:] < C[1:, i0 - 1:].min(0))

# 1605: last digits 1 and 9 lead 3 and 7 only 1 to 3 % of the time
Q = P[P > 5]; C = {c: np.cumsum(Q % 10 == c) for c in (1, 3, 7, 9)}
assert [round(100 * float(np.mean(C[a] > C[b])), 1) for a, b in ((9, 3), (9, 7), (1, 3), (1, 7))] == [1.0, 1.2, 2.5, 2.8]

# 1606: first leads of the Kronecker races for k = -4, 12, 3, 2
def kron(k, q):
    if q == 2: return 0 if k % 2 == 0 else (1 if k % 8 in (1, 7) else -1)
    return legendre_symbol(k % q, q) if k % q else 0
def first_lead(k, X):
    Pk = P[P < X]; per = 4 * abs(k); small = Pk[Pk <= per]; big = Pk[Pk > per]
    tab = np.zeros(per, np.int64)
    for r in range(per):
        if math.gcd(r, per) == 1: tab[r] = kron(k, int(big[np.argmax(big % per == r)]))
    c = np.cumsum(np.concatenate(([kron(k, int(q)) for q in small], tab[big % per])))
    return int(Pk[np.argmax(c > 0)])
assert [first_lead(k, 1.2e7) for k in (-4, 12, 3, 2)] == [26861, 61463, 61981, 11100143]

# 1606: ties in the mod 4 race
c4 = np.cumsum(P % 4 == 1) - np.cumsum(P % 4 == 3)
assert [int(x) for x in P[c4 == 0][:6]] == [2, 5, 17, 41, 461, 26833]

# 1305: record least quadratic non-residues below 10^6
def records(f, X):
    out, best = [], 0
    for p in primerange(3, X):
        v = f(p)
        if v > best: best = v; out.append((v, p))
    return out
def n0(p):
    a = 2
    while legendre_symbol(a, p) == 1: a += 1
    return a
assert records(n0, 10**6) == [(2, 3), (3, 7), (5, 23), (7, 71), (11, 311), (13, 479), (17, 1559), (19, 5711), (23, 10559), (29, 18191), (31, 31391), (43, 366791)]

# 1404: record least primitive roots below 10^6
def g0(p):
    qs = primefactors(p - 1); a = 2
    while any(pow(a, (p - 1) // q, p) == 1 for q in qs): a += 1
    return a
assert records(g0, 10**6) == [(2, 3), (3, 7), (5, 23), (6, 41), (7, 71), (19, 191), (21, 409), (23, 2161), (31, 5881), (37, 36721), (38, 55441), (44, 71761), (69, 110881), (73, 760321)]

# 1211: Wieferich primes to the bases 3, 5, 7, 10, 11
P7 = [int(p) for p in P[P < 10**7]]
for b, known in ((3, [11, 1006003]), (5, [2, 20771, 40487, 53471161]), (7, [5, 491531]), (10, [3, 487, 56598313]), (11, [71])):
    assert [p for p in P7 if b % p and gmpy2.powmod(b, p - 1, p * p) == 1] == [x for x in known if x < 10**7]
    assert all(pow(b, x - 1, x * x) == 1 and isprime(x) for x in known if x > 2)

# 1212: Giuga's criterion for 2 <= n < 2000
assert [n for n in range(2, 2000) if sum(pow(i, n - 1, n) for i in range(1, n)) % n == n - 1] == list(primerange(2, 2000))

# 1509: Zsigmondy exceptions for a^n - 1 with a <= 30, n <= 60
exc = set()
for a in range(2, 31):
    for n in range(2, 61):
        v = int(cyclotomic_poly(n, a))
        for q in primefactors(n):
            while v % q == 0: v //= q
        if v == 1: exc.add((a, n))
assert exc == {(2, 6), (3, 2), (7, 2), (15, 2)}

# 1511: Catalan, 8 and 9 are the only consecutive perfect powers below 10^12
pw = sorted({b ** e for e in range(2, 41) for b in range(2, int(round(10 ** (12 / e))) + 2) if b ** e < 10**12})
assert [x for x, y in zip(pw, pw[1:]) if y - x == 1] == [8]

# 1513: Sylvester's sequence, pairwise coprime, sum of reciprocals, no factor = 5 (mod 6)
sy = [2]
for _ in range(7): sy.append(sy[-1] ** 2 - sy[-1] + 1)
assert sy[:6] == [2, 3, 7, 43, 1807, 3263443] and all(math.gcd(a, b) == 1 for i, a in enumerate(sy) for b in sy[:i])
assert all(Fraction(1, 1) - sum(Fraction(1, x) for x in sy[:n]) == Fraction(1, sy[n] - 1) for n in range(1, 8))
assert all(q % 6 != 5 for x in sy[:7] for q in factorint(x))

# 1609: Chebotarev for x^3 - 2, the shares of 3, 1, 0 roots below 10^6
cnt = [0, 0, 0]
for p in primerange(5, 10**6):
    cnt[1 if p % 3 == 2 else (0 if pow(2, (p - 1) // 3, p) == 1 else 2)] += 1
assert [round(c / sum(cnt), 4) for c in cnt] == [0.1660, 0.5002, 0.3338]

# 1607: the least prime in each class mod q, against (phi(q) ln q)^2, 3 < q <= 300
lin = 0
for q in range(4, 301):
    r = P[:200000] % q; u, idx = np.unique(r, return_index=True); ok = np.gcd(u, q) == 1
    assert ok.sum() == totient(q)
    lin = max(lin, int(P[idx[ok]].max()) / (int(totient(q)) * math.log(q)) ** 2)
assert 0.650 < lin < 0.651

# 1510: the least prime = 1 (mod d) is at most 2^d - 1
assert all(2 ** d - 1 >= min(int(p) for p in P[:50000] if p % d == 1) for d in range(3, 200))

# 1702, 1703: irregular primes below 2000, via p B_k = sum j^k (mod p^2)
irr = {}
for p in primerange(5, 2000):
    m = p * p; j = np.arange(1, p, dtype=np.int64); j2 = j * j % m; pwr = j2.copy(); ks = []
    for k in range(2, p - 2, 2):
        if int(pwr.sum() % m) == 0: ks.append(k)
        pwr = pwr * j2 % m
    if ks: irr[p] = ks
assert list(irr)[:12] == [37, 59, 67, 101, 103, 131, 149, 157, 233, 257, 263, 271] and sum(p < 1000 for p in irr) == 64
assert len(irr) == 121 and len(list(primerange(5, 2000))) == 301 and irr[157] == [62, 110] and irr[491] == [292, 336, 338]
assert min(p for p in irr if len(irr[p]) == 3) == 491 and min(p for p in irr if len(irr[p]) == 2) == 157

# 1704: 16843 is a Wolstenholme prime, C(2p-1, p-1) = 1 (mod p^4)
p = 16843; m = p**4; c = 1
for i in range(1, p): c = c * (p + i) % m * pow(i, -1, m) % m     # C(2p-1, p-1) = prod (p+i)/i
assert c == 1

# 2006, 2008, 2024: Mertens and Liouville sums to 10^7, and the mean of omega(n)
X = 10**7; mu = np.ones(X + 1, np.int8); lam = np.ones(X + 1, np.int8); om = np.zeros(X + 1, np.int8)
for p in P[P <= X].tolist():
    mu[::p] *= -1; om[::p] += 1
    if p * p <= X: mu[::p * p] = 0
    pk = p
    while pk <= X: lam[::pk] *= -1; pk *= p
mu[0] = lam[0] = 0
M = np.cumsum(mu, dtype=np.int64); Lv = np.cumsum(lam, dtype=np.int64); sq = np.sqrt(np.arange(X + 1))
assert [int(M[10**k]) for k in range(1, 8)] == [-1, 1, 2, -23, -48, 212, 1037] and np.all(np.abs(M[2:]) < sq[2:])
rat = np.abs(M[101:]) / sq[101:]; assert int(np.argmax(rat)) + 101 == 199 and round(float(rat.max()), 3) == 0.567
assert [int(Lv[10**k]) for k in range(1, 8)] == [0, -2, -14, -94, -288, -530, -842] and Lv[2:].max() == 0
assert round(float(om[1:10**6 + 1].mean()), 3) == 2.854
del mu, lam, om, M, Lv, sq

# 2009: Robin's inequality fails for exactly 26 values n >= 3, the last 5040, and holds on (5040, 10^6]
Y = 10**6; sig = np.zeros(Y + 1, np.int64)
for d in range(1, Y + 1): sig[d::d] += d
nn = np.arange(3, Y + 1); viol = nn[sig[3:] >= np.exp(0.5772156649015329) * nn * np.log(np.log(nn))]
assert viol.tolist() == [3, 4, 5, 6, 8, 9, 10, 12, 16, 18, 20, 24, 30, 36, 48, 60, 72, 84, 120, 180, 240, 360, 720, 840, 2520, 5040]
assert round(19344 / (math.exp(0.5772156649015329) * 5040 * math.log(math.log(5040))), 4) == 1.0056 and sig[5040] == 19344

# 2009: Lagarias' inequality for 1 < n <= 10^5
H = np.cumsum(1 / np.arange(1, 10**5 + 1))
assert np.all(sig[2:10**5 + 1] < H[1:] + np.exp(H[1:]) * np.log(H[1:]))

# 6107: no composite n <= 10^6 has phi(n) | n - 1
phi = np.arange(Y + 1, dtype=np.int64)
for p in P[P <= Y].tolist(): phi[::p] -= phi[::p] // p
comp = np.nonzero(~isp[:Y + 1])[0]; comp = comp[comp > 1]
assert not np.any((comp - 1) % phi[comp] == 0)

# 2022: records of the least prime in a Goldbach partition, n <= 1.1e7
Z = 11_000_000; ev = np.arange(6, Z + 1, 2, dtype=np.int64); least = np.zeros(len(ev), np.int32); op = np.arange(len(ev))
for p in P[1:].tolist():
    hit = isp[ev[op] - p]; least[op[hit]] = p; op = op[~hit]
    if not len(op): break
rec, best = [], 0
for i in np.nonzero(least > np.maximum.accumulate(np.concatenate(([0], least[:-1]))))[0]:
    if least[i] > best: best = int(least[i]); rec.append((int(ev[i]), best))
assert rec == [(6, 3), (12, 5), (30, 7), (98, 19), (220, 23), (308, 31), (556, 47), (992, 73), (2642, 103), (5372, 139), (7426, 173),
               (43532, 211), (54244, 233), (63274, 293), (113672, 313), (128168, 331), (194428, 359), (194470, 383), (413572, 389),
               (503222, 523), (1077422, 601), (3526958, 727), (3807404, 751), (10759922, 829)]

# 2022: the 2012 record partition 3325581707333960528 = 9781 + p, and no smaller prime works
gn = 3325581707333960528
assert isprime(gn - 9781) and not any(isprime(q) and isprime(gn - q) for q in range(3, 9781, 2))

# 5007: the number of n <= 10^k with n^2 + 1 prime
c2, out, nxt = 0, [], 10
for n in range(1, 10**6 + 1):
    c2 += gmpy2.is_prime(n * n + 1)
    if n == nxt: out.append(c2); nxt *= 10
assert out == [5, 19, 112, 841, 6656, 54110]

# 6104: base-2 Fermat and strong pseudoprimes and Carmichael numbers below 10^k
def powmod_v(b, e, m):
    r = np.ones_like(m); b = b % m; e = e.copy()
    while e.any():
        r = np.where(e & 1 == 1, r * b % m, r); b = b * b % m; e >>= 1
    return r
oc = np.arange(3, 10**7, 2, dtype=np.int64); oc = oc[~isp[oc]]
psp = oc[powmod_v(np.full_like(oc, 2), oc - 1, oc) == 1]
assert [int((psp <= 10**k).sum()) for k in range(3, 8)] == [3, 22, 78, 245, 750]
sp = [int(n) for n in psp if gmpy2.is_strong_prp(int(n), 2)]
assert [sum(x <= 10**k for x in sp) for k in range(3, 8)] == [0, 5, 16, 46, 162] and sp[0] == 2047
car = [int(n) for n in psp if all(e == 1 for e in factorint(int(n)).values()) and all((n - 1) % (q - 1) == 0 for q in factorint(int(n)))]
assert [sum(x <= 10**k for x in car) for k in range(3, 8)] == [1, 7, 16, 43, 105]

# 6105: the least strong pseudoprimes to the bases 2, 3 and 2, 3, 5, and the table for the first 13 prime bases
assert min(n for n in sp if gmpy2.is_strong_prp(n, 3)) == 1373653 and not any(gmpy2.is_strong_prp(n, 3) and gmpy2.is_strong_prp(n, 5) for n in sp)
mr = [2047, 1373653, 25326001, 3215031751, 2152302898747, 3474749660383, 341550071728321, 341550071728321,
      3825123056546413051, 3825123056546413051, 3825123056546413051, 318665857834031151167461, 3317044064679887385961981]
for n, x in enumerate(mr, 1):
    assert not isprime(x) and all(gmpy2.is_strong_prp(x, b) for b in list(primerange(2, 50))[:n])
    assert n == 13 or not all(gmpy2.is_strong_prp(x, b) for b in list(primerange(2, 50))[:n + 1]) or mr[n] == x

# 6104: Chernick's (6k+1)(12k+1)(18k+1) is a Carmichael number when the factors are prime
for k in range(1, 100):
    a, b, c = 6 * k + 1, 12 * k + 1, 18 * k + 1
    if isprime(a) and isprime(b) and isprime(c): assert all((a * b * c - 1) % (q - 1) == 0 for q in (a, b, c))

# 6202: Euler-Lagrange, for p = 3 (mod 4), 2p + 1 is prime iff it divides 2^p - 1
assert all((pow(2, p, 2 * p + 1) == 1) == isprime(2 * p + 1) for p in primerange(5, 10**5) if p % 4 == 3)

# 6203: covering sets of the Riesel number 509203 and the prime Sierpinski number 271129
for k, cover in ((509203, (3, 5, 7, 13, 17, 241)), (271129, (3, 5, 7, 13, 17, 241))):
    sgn = -1 if k == 509203 else 1
    assert all(any((k * pow(2, n, q) + sgn) % q == 0 for q in cover) for n in range(1, 200))

# 6204: factorial primes n! +- 1 for n <= 500
f, fp, fm = 1, [], []
for n in range(1, 501):
    f *= n
    if gmpy2.is_prime(f + 1): fp.append(n)
    if gmpy2.is_prime(f - 1): fm.append(n)
assert fp == [1, 2, 3, 11, 27, 37, 41, 73, 77, 116, 154, 320, 340, 399, 427]
assert fm == [3, 4, 6, 7, 12, 14, 30, 32, 33, 38, 94, 166, 324, 379, 469]

# 6205: primorial primes p# +- 1 for p < 2000
pr, pp, pm = 1, [], []
for p in primerange(2, 2000):
    pr *= p
    if gmpy2.is_prime(pr + 1): pp.append(p)
    if gmpy2.is_prime(pr - 1): pm.append(p)
assert pp == [2, 3, 5, 7, 11, 31, 379, 1019, 1021] and pm == [3, 5, 11, 13, 41, 89, 317, 337, 991, 1873]

# 6206: the first 100 Fortunate numbers are prime
pr, fo = 1, []
for p in primerange(2, 600):
    pr *= p; fo.append(int(gmpy2.next_prime(pr + 1)) - pr)
assert fo[:10] == [3, 5, 7, 13, 23, 17, 19, 23, 37, 61] and all(isprime(x) for x in fo[:100])

# 6204: Brocard, n! + 1 is a square only for n = 4, 5, 7 (n <= 1000)
f = 1; assert [n for n in range(1, 1001) if gmpy2.is_square((f := f * n) + 1)] == [4, 5, 7]

# 6205: the first 9 terms of the Euclid-Mullin sequence
em = [2]
while len(em) < 9: em.append(min(factorint(math.prod(em) + 1)))
assert em == [2, 3, 7, 43, 13, 53, 5, 6221671, 38709183810571]

# 6207: the 46 Pierpont primes below 2e6
pier = sorted(2**u * 3**v + 1 for u in range(22) for v in range(14) if 2**u * 3**v + 1 < 2 * 10**6 and isprime(2**u * 3**v + 1))
assert len(pier) == 46 and pier[:12] == [2, 3, 5, 7, 13, 17, 19, 37, 73, 97, 109, 163]

# 6209: Wagstaff exponents below 1000
assert [p for p in primerange(3, 1000) if gmpy2.is_prime((2**p + 1) // 3)] == [3, 5, 7, 11, 13, 17, 19, 23, 31, 43, 61, 79, 101, 127, 167, 191, 199, 313, 347, 701]

# 6210: repunit primes R_n for n < 1100
assert [n for n in range(2, 1100) if gmpy2.is_prime((10**n - 1) // 9)] == [2, 19, 23, 317, 1031]

# 6208: Cullen and Woodall primes for n <= 1000
assert [n for n in range(1, 1001) if gmpy2.is_prime(n * 2**n + 1)] == [1, 141] and [n for n in range(1, 1001) if gmpy2.is_prime(n * 2**n - 1)] == [2, 3, 6, 30, 75, 81, 115, 123, 249, 362, 384, 462, 512, 751, 822]

# 6211: palindromic primes with d <= 9 digits
def pal(d):
    h = (d + 1) // 2; out = 0
    for x in range(10**(h - 1), 10**h):
        s_ = str(x); v = int(s_ + s_[::-1][d % 2:])
        out += bool(isp[v]) if v <= 200_000_000 else isprime(v)
    return out
assert [pal(d) for d in range(1, 10)] == [4, 1, 15, 0, 93, 0, 668, 0, 5172]

# 6212: reversible primes with d <= 7 digits
def rev(a):
    r = np.zeros_like(a); a = a.copy()
    while a.any(): r = r * 10 + a % 10; a //= 10
    return r
assert [int(at(rev(P[(P >= 10**(d - 1)) & (P < 10**d)])).sum()) for d in range(1, 8)] == [4, 9, 43, 204, 1499, 9538, 71142]

# 6211: Belphegor-type primes 1 0^k 666 0^k 1 for k <= 100
assert [k for k in range(101) if gmpy2.is_prime(10**(2 * k + 4) + 666 * 10**(k + 1) + 1)] == [0, 13, 42]

# 6213: Gaussian primes with |z| <= 10 and |z| <= 100
gp = lambda R: sum(1 for a in range(-R, R + 1) for b in range(-R, R + 1) if a * a + b * b <= R * R and
                   (isprime(a * a + b * b) if a and b else (isprime(abs(a + b)) and abs(a + b) % 4 == 3)))
assert gp(10) == 100 and gp(100) == 4928

# 6303: the least titanic prime is 10^999 + 7
assert int(gmpy2.next_prime(10**999)) == 10**999 + 7

# --- second research pass (October 2026) --------------------------------------------------------------

# 3019: first occurrence of each gap below 2e8
u, idx = np.unique(g, return_index=True); fo = {int(a): int(P[i]) for a, i in zip(u, idx)}
assert [fo[k] for k in range(2, 42, 2)] == [3, 7, 23, 89, 139, 199, 113, 1831, 523, 887, 1129, 1669, 2477, 2971, 4297, 5591, 1327, 9551, 30593, 19333]
assert all(k in fo for k in range(2, 200, 2)) and 200 not in fo and max(fo) == 248

# 3020: primes between n^2 and (n+1)^2, Oppermann, Brocard
pi_ = lambda x: np.searchsorted(P, np.asarray(x), side='right')
nn = np.arange(1, 14142); leg = pi_((nn + 1) ** 2) - pi_(nn ** 2)
assert leg.min() == 2 and [int(nn[leg <= k].max()) for k in range(2, 8)] == [5, 9, 14, 17, 23, 26]
lo, hi = pi_(nn[1:] ** 2 - 1) - pi_(nn[1:] ** 2 - nn[1:]), pi_(nn[1:] ** 2 + nn[1:]) - pi_(nn[1:] ** 2)
assert lo.min() >= 1 and hi.min() >= 1 and int(nn[1:][lo == 1].max()) == 11 and int(nn[1:][hi == 1].max()) == 17
bp = P[(P >= 3) & (P < 14142)]
assert int((pi_(bp[1:] ** 2) - pi_(bp[:-1] ** 2)).min()) == 5

# 1610: every class mod q has a prime below q^2, and max P(q)/(q ln^2 q) is at q = 461 (q <= 500)
lmax, lq = 0, 0
for q in range(3, 501):
    r = P[:300000] % q; uu, ii = np.unique(r, return_index=True); ok = np.gcd(uu, q) == 1
    assert ok.sum() == totient(q)
    far = int(P[ii[ok]].max()); assert far < q * q
    if far / (q * math.log(q) ** 2) > lmax: lmax, lq = far / (q * math.log(q) ** 2), (q, far)
assert lq == (461, 37363) and round(lmax, 3) == 2.154

# 1213: no Wall-Sun-Sun prime below 10^6
def fib_mod(n, m):                                 # (F_n, F_(n+1)) mod m by fast doubling
    if n == 0: return 0, 1
    a, b = fib_mod(n // 2, m); c = a * (2 * b - a) % m; d = (a * a + b * b) % m
    return (d, (c + d) % m) if n % 2 else (c, d)
assert not any(fib_mod(p - legendre_symbol(5, p), p * p)[0] == 0 for p in primerange(7, 10**6))

# 1514: Fibonacci, Lucas and Pell primes for n <= 1000
assert [n for n in range(1, 1001) if gmpy2.is_prime(int(fibonacci(n)))] == [3, 4, 5, 7, 11, 13, 17, 23, 29, 43, 47, 83, 131, 137, 359, 431, 433, 449, 509, 569, 571]
assert [n for n in range(0, 1001) if gmpy2.is_prime(int(lucas(n)))] == [0, 2, 4, 5, 7, 8, 11, 13, 16, 17, 19, 31, 37, 41, 47, 53, 61, 71, 79, 113, 313, 353, 503, 613, 617, 863]
pell = [0, 1]
for _ in range(1000): pell.append(2 * pell[-1] + pell[-2])
assert [n for n in range(1001) if gmpy2.is_prime(pell[n])] == [2, 3, 5, 11, 13, 29, 41, 53, 59, 89, 97, 101, 167, 181, 191, 523, 929]

# 1514: F_(2n+1) +- 2 factor as F_(n+2) L_(n-1) and F_(n-1) L_(n+2)
assert all(sorted([fibonacci(2 * n + 1) + 2, fibonacci(2 * n + 1) - 2]) == sorted([fibonacci(n + 2) * lucas(n - 1), fibonacci(n - 1) * lucas(n + 2)]) for n in range(2, 300))

# 1515: de Polignac numbers below 4000 and Pollack's covering congruence
pol = [n for n in range(1, 4000, 2) if not any(isprime(n - 2**k) for k in range(n.bit_length()) if n - 2**k > 1)]
assert pol[:10] == [1, 127, 149, 251, 331, 337, 373, 509, 599, 701] and min(n for n in pol if not isprime(n) and n > 1) == 905 and 959 in pol
assert all(any((1260327937 - pow(2, k, q)) % q == 0 for q in (3, 5, 7, 13, 17, 241)) for k in range(24)) and 2863311360 == 2**9 * 3 * 5 * 7 * 13 * 17 * 241

# 1516: 2^k - 7 is prime for k = 39, 715, 1983 (k <= 2000)
assert [k for k in range(3, 2001) if gmpy2.is_prime(2**k - 7)] == [39, 715, 1983]

# 1518: Reyssat's abc triple has quality 1.62991
a_, b_, c_ = 2, 3**10 * 109, 23**5; assert a_ + b_ == c_ and round(math.log(c_) / math.log(2 * 3 * 109 * 23), 5) == 1.62991

# 1114: 123...n is coprime to 30 only for n = 1, 7, 13, 19 (mod 30)
assert all(n % 30 in (1, 7, 13, 19) for n in range(1, 400) if int(''.join(map(str, range(1, n + 1)))) % 2 and int(''.join(map(str, range(1, n + 1)))) % 3 and int(''.join(map(str, range(1, n + 1)))) % 5)

# 2033: Dusart's bounds with sharp thresholds 599, 88789 (not 88783), and the 1.2762 upper bound
lnP = np.log(P.astype(np.float64)); xs = P.astype(np.float64); idx = np.arange(1, len(P) + 1, dtype=np.float64)
f1 = lambda x: x / np.log(x) * (1 + 1 / np.log(x))                 # holds for x >= 599, fails just below 599
assert P[np.nonzero(idx < f1(xs))[0][-1]] < 599 and P[np.nonzero(idx - 1 < f1(xs - 1e-9))[0][-1]] == 599
fb = lambda x: x / np.log(x) * (1 + 1 / np.log(x) + 2 / np.log(x) ** 2)
assert int(np.searchsorted(P, 88783, 'right')) < fb(88783) and all(int(np.searchsorted(P, x, 'right')) >= fb(x) for x in range(88789, 200000))
assert np.all(idx[1:] <= xs[1:] / lnP[1:] * (1 + 1.2762 / lnP[1:]))

# 2034: theta(x) < x below 2e8, max theta(x)/x = 0.999976
th = np.cumsum(lnP)
assert np.all(th < xs) and round(float((th / xs)[100:].max()), 6) == 0.999976

# 2033: |theta(x) - x| <= 0.2 x / ln^2 x from 3594641 on, failing just below it
i0 = int(np.searchsorted(P, 3594641)); bnd = 0.2 * xs / lnP ** 2
assert np.all(np.abs(th[i0:] - xs[i0:]) <= bnd[i0:]) and np.all(xs[i0 + 1:] - th[i0:-1] <= bnd[i0 + 1:])
assert xs[i0] - th[i0 - 1] > bnd[i0] and P[i0 - 1] == 3594631

# 2033: a prime in [x, x + x/(25 ln^2 x)] for x >= 396738, the last failure is the gap after 396733
h25 = xs[:-1] / (25 * lnP[:-1] ** 2); j0 = int(np.searchsorted(P, 396738))
assert np.all(g[j0:] <= h25[j0:]) and P[j0] <= 396738 * (1 + 1 / (25 * math.log(396738) ** 2)) and int(P[np.nonzero(g > h25)[0][-1]]) == 396733

# 2035: the Landau-Ramanujan constant from the product over p = 3 (mod 4), p < 10^7
Q3 = P[(P % 4 == 3) & (P < 10**7)].astype(np.float64)
assert abs(float(np.prod((1 - Q3 ** -2) ** -0.5)) / math.sqrt(2) - 0.7642236535) < 3e-9

# 3021: extreme ratios of consecutive gaps below 2e8 (104 and 1/104)
r_ = g[1:] / g[:-1]; assert r_.max() == 104 and int(P[np.argmax(r_)]) == 192983849 and int(P[np.argmin(r_)]) == 195208231

# 4013: the gaps 2, 4, 8, 16, 32 first follow each other at 6824897
w = [2, 4, 8, 16, 32]; ok = np.ones(len(g) - 5, bool)
for j, v in enumerate(w): ok &= g[j:len(g) - 5 + j] == v
assert int(P[np.argmax(ok)]) == 6824897

# 4014: balanced, strong and weak primes; the n-th strong prime is first below the n-th weak one at n = 19799
mid = P[1:-1]; bal = mid[2 * mid == P[:-2] + P[2:]]; assert bal[:10].tolist() == [5, 53, 157, 173, 211, 257, 263, 373, 563, 593]
stg = mid[2 * mid > P[:-2] + P[2:]]; wk = mid[2 * mid < P[:-2] + P[2:]]
assert stg[:4].tolist() == [11, 17, 29, 37] and wk[:5].tolist() == [3, 7, 13, 19, 23] and int(np.argmax(stg[:30000] < wk[:30000])) + 1 == 19799

# 4014: the first good primes
gp_ = [int(P[n]) for n in range(1, 40) if all(int(P[n]) ** 2 > int(P[n - i]) * int(P[n + i]) for i in range(1, n + 1))]
assert gp_[:10] == [5, 11, 17, 29, 37, 41, 53, 59, 67, 71]

# 4015: Chen primes below 100
chen = [p for p in primerange(2, 100) if sum(factorint(p + 2).values()) <= 2]
assert chen[:14] == [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 47] and [p for p in primerange(2, 100) if p not in chen] == [43, 61, 73, 79, 97]

# 4016: every integer above 96 (up to 20000) is a sum of distinct super-primes
sup = [int(P[int(p) - 1]) for p in P[:3000]]; reach = 1
for s_ in sup: reach |= reach << s_
reach &= (1 << 20001) - 1
assert sup[:10] == [3, 5, 11, 17, 31, 41, 59, 67, 83, 109] and max(n for n in range(20001) if not reach >> n & 1) == 96

# 4016: the first lucky primes
lk = list(range(1, 2000, 2)); i = 1
while i < len(lk) and lk[i] <= len(lk):
    st = lk[i]; del lk[st - 1::st]; i += 1
assert [x for x in lk if isprime(x)][:13] == [3, 7, 13, 31, 37, 43, 67, 73, 79, 127, 151, 163, 193]

# 6221: Pillai primes below 80
pil = [p for p in primerange(2, 80) if any((math.factorial(m) + 1) % p == 0 and p % m != 1 for m in range(2, p))]
assert pil == [23, 29, 59, 61, 67, 71, 79]

# 6221: self primes below 250
dsum = lambda n: sum(map(int, str(n)))
gen = {n + dsum(n) for n in range(1, 3000)}; assert [p for p in primerange(2, 250) if p not in gen] == [3, 5, 7, 31, 53, 97, 211, 233]

# 6221: primeval numbers, 1379 gives 31 primes and 13679 gives 106
def primeval(n):
    d = str(n); found = set()
    for r in range(1, len(d) + 1):
        for t in permutations(d, r):
            if t[0] != '0' and isprime(int(''.join(t))): found.add(int(''.join(t)))
    return len(found)
assert primeval(1379) == 31 and primeval(13679) == 106

# 6215: truncatable primes, 4260 left, 83 right, 15 both
def trunc(left):
    layer, tot = [2, 3, 5, 7], [2, 3, 5, 7]
    while layer:
        layer = [int(d + str(p)) if left else p * 10 + int(d) for p in layer for d in ('123456789' if left else '1379')]
        layer = [q for q in layer if isprime(q)]; tot += layer
    return tot
lt, rt = trunc(True), trunc(False)
assert len(lt) == 4260 and max(lt) == 357686312646216567629137 and len(rt) == 83 and max(rt) == 73939133
assert sorted(set(lt) & set(rt)) == [2, 3, 5, 7, 23, 37, 53, 73, 313, 317, 373, 797, 3137, 3797, 739397]

# 6216: circular and permutable primes below 10^7
P7s = set(int(p) for p in P[P < 10**7])
circ = sorted({min(int(str(p)[i:] + str(p)[:i]) for i in range(len(str(p)))) for p in P7s
               if all(int(str(p)[i:] + str(p)[:i]) in P7s for i in range(len(str(p))))})
assert circ == [2, 3, 5, 7, 11, 13, 17, 37, 79, 113, 197, 199, 337, 1193, 3779, 11939, 19937, 193939, 199933]
perm = sorted(p for p in P7s if set(str(p)) <= set('1379') or p < 10 if all(int(''.join(t)) in P7s for t in set(permutations(str(p)))))
assert perm == [2, 3, 5, 7, 11, 13, 17, 31, 37, 71, 73, 79, 97, 113, 131, 199, 311, 337, 373, 733, 919, 991]

# 6217: the 26 minimal primes are prime, pairwise incomparable, and cover every prime below 10^6
mins = ['2', '3', '5', '7', '11', '19', '41', '61', '89', '409', '449', '499', '881', '991', '6469', '6949', '9001', '9049', '9649', '9949', '60649', '666649', '946669', '60000049', '66000049', '66600049']
def subseq(a, b):
    it = iter(b); return all(c in it for c in a)
assert all(isprime(int(m)) for m in mins) and not any(subseq(a, b) for a in mins for b in mins if a != b)
assert all(any(subseq(m, str(p)) for m in mins) for p in P[P < 10**6].tolist())

# 6218: happy primes below 10^4, 10^5, 10^6, and M42643801 is happy
def happy(n):
    seen = set()
    while n != 1 and n not in seen: seen.add(n); n = sum(int(c) ** 2 for c in str(n))
    return n == 1
hp = [int(p) for p in P[P < 10**6] if happy(int(p))]; assert [sum(x < 10**k for x in hp) for k in (4, 5, 6)] == [200, 1465, 11144]
assert happy(sum(int(c) ** 2 for c in (gmpy2.mpz(2) ** 42643801 - 1).digits()))

# 6218: 3294 Smith numbers below 10^5
spf = np.zeros(10**5 + 1, np.int64)
for p in P[P <= 316].tolist():
    m_ = spf[p * p::p] == 0; spf[p * p::p][m_] = p
def factors(n):
    out = []
    while spf[n]: out.append(int(spf[n])); n //= int(spf[n])
    return out + [n]
assert sum(1 for n in range(4, 10**5) if spf[n] and dsum(n) == sum(dsum(q) for q in factors(n))) == 3294

# 6218: the home prime HP(8)
def home(n):
    while not isprime(n): n = int(''.join(str(q) * e for q, e in sorted(factorint(n).items())))
    return n
assert home(8) == 3331113965338635107

# 6219: strobogrammatic primes below 70000
strob = {'0': '0', '1': '1', '8': '8', '6': '9', '9': '6'}
assert [p for p in P[P < 70000].tolist() if set(str(p)) <= set(strob) and int(''.join(strob[c] for c in str(p)[::-1])) == p] == [11, 101, 181, 619, 16091, 18181, 19861, 61819]

# 6219: the least and largest primes using the digits 2, 3, 5, 7 once
pd = [p for p in P[P < 10000].tolist() if sorted(str(p)) == ['2', '3', '5', '7']]; assert pd[0] == 2357 and pd[-1] == 7523

# 6219: decimal covering sets for 91...13, 9175 * 10^n + 1 and 10176 * 10^n - 1
assert all(any(((82 * 10**n + 17) // 9) % q == 0 for q in (3, 7, 11, 13)) for n in range(1, 200))
assert all(any((9175 * 10**n + 1) % q == 0 for q in (7, 11, 13, 37)) and any((10176 * 10**n - 1) % q == 0 for q in (7, 11, 13, 37)) for n in range(1, 200))

# 6220: prime prefixes of pi and e up to 1000 digits
mpmath.mp.dps = 1020
for const, lens in ((mpmath.pi, [1, 2, 6, 38]), (mpmath.e, [1, 3, 7, 85])):
    dig = mpmath.nstr(const, 1010, strip_zeros=False).replace('.', '')
    assert [k for k in range(1, 1001) if gmpy2.is_prime(int(dig[:k]))] == lens

# 6220: the first 10-digit prime in e is 7427466391, at the 99th decimal
dig = mpmath.nstr(mpmath.e, 1010).replace('.', '')
assert next(i for i in range(1000) if gmpy2.is_prime(int(dig[i:i + 10])) and dig[i] != '0') == 99 and dig[99:109] == '7427466391'

# 6220: the continued fraction of the Copeland-Erdos constant starts [0; 4, 4, 8, 16, 18, 5, 1, 1, 1, 1, 7]
ce = ''.join(str(p) for p in primerange(2, 2000)); lo_ = Fraction(int(ce[:300]), 10**300)
cfr, x_ = [], lo_
for _ in range(12):
    a = x_.numerator // x_.denominator; cfr.append(a); x_ = 1 / (x_ - a)
assert cfr == [0, 4, 4, 8, 16, 18, 5, 1, 1, 1, 1, 7]

# 2029: N(100) = 29 and N(1000) = 649 zeros of zeta
mpmath.mp.dps = 20
assert mpmath.nzeros(100) == 29 and mpmath.nzeros(1000) == 649

# --- third research pass (October 2026) ---------------------------------------------------------------

# 4017: sums of the first n primes that are prime, divisible by n, or a square (n <= 4e6)
S_ = np.cumsum(P[:4_000_000].astype(np.int64)); nn_ = np.arange(1, len(S_) + 1)
assert [int(n) for n in nn_[:200] if isprime(int(S_[n - 1]))][:10] == [1, 2, 4, 6, 12, 14, 60, 64, 96, 100]
assert nn_[S_ % nn_ == 0].tolist() == [1, 23, 53, 853, 11869, 117267, 339615, 3600489]
rt = np.round(np.sqrt(S_.astype(np.float64))).astype(np.int64)
assert nn_[rt * rt == S_].tolist() == [9, 2474, 6694, 7785, 709838] and int(P[709837]) == 10729219 and int(rt[709837]) == 1916357

# 4017: 3 + 5 + ... + p_k is a square for k = 24, 28, 32, 46, 296, 5327, 12227 (k = index of the last prime)
O_ = np.cumsum(P[1:].astype(np.int64)); ro = np.round(np.sqrt(O_.astype(np.float64))).astype(np.int64)
assert (np.nonzero(ro * ro == O_)[0] + 2).tolist() == [24, 28, 32, 46, 296, 5327, 12227]

# 4018: the least numbers with exactly m = 1..8 representations as sums of consecutive primes
Xc = 4_000_000; Pc = P[P <= Xc].astype(np.int64); C0 = np.concatenate(([0], np.cumsum(Pc))); ways = np.zeros(Xc + 1, np.int32)
for Lk in range(1, len(Pc) + 1):
    sm = C0[Lk:] - C0[:-Lk]; sm = sm[sm <= Xc]
    if not len(sm): break
    np.add.at(ways, sm, 1)
assert [int(np.argmax(ways == m)) for m in range(1, 9)] == [2, 5, 41, 1151, 311, 34421, 218918, 3634531]

# 4019: least prime k-tuplets of minimal width for k = 5..8 (the k primes consecutive)
def pattern_starts(h, X=2 * 10**8):
    ok = np.ones(len(P), bool)
    for d in h[1:]: ok &= at(P + d)
    st = P[ok]
    return [int(p) for p in st if all(not at(np.array([p + x]))[0] for x in range(1, h[-1]) if x not in h)]
tup_tab = {(0, 2, 6, 8, 12): [5, 11, 101], (0, 4, 6, 10, 12): [7, 97, 1867], (0, 4, 6, 10, 12, 16): [7, 97, 16057],
           (0, 2, 6, 8, 12, 18, 20): [11, 165701, 1068701], (0, 2, 8, 12, 14, 18, 20): [5639, 88799, 284729],
           (0, 2, 6, 8, 12, 18, 20, 26): [11, 15760091, 25658441], (0, 2, 6, 12, 14, 20, 24, 26): [17, 1277, 113147],
           (0, 6, 8, 14, 18, 20, 24, 26): [88793, 284723, 855713]}
for h, first3 in tup_tab.items(): assert pattern_starts(h)[:3] == first3

# 4019: quadruplets are 30n + 11.. with centre 15, 105, 195 (mod 210); counts by digits
q4 = P[at(P + 2) & at(P + 6) & at(P + 8)]
assert all(int(p) % 30 == 11 and (int(p) + 4) % 210 in (15, 105, 195) for p in q4[1:])
assert [int(((q4 >= 10**(d - 1)) & (q4 < 10**d)).sum()) for d in range(2, 9)] == [1, 3, 7, 26, 128, 733, 3869]

# 4019: quintuplet and sextuplet counts below 10^k
qa, qb = P[at(P + 2) & at(P + 6) & at(P + 8) & at(P + 12)], P[at(P + 4) & at(P + 6) & at(P + 10) & at(P + 12)]
sx = P[at(P + 4) & at(P + 6) & at(P + 10) & at(P + 12) & at(P + 16)]
assert [int((qa <= 10**k).sum()) for k in range(1, 9)] == [1, 2, 3, 4, 10, 34, 160, 697]
assert [int((qb <= 10**k).sum()) for k in range(1, 9)] == [1, 2, 2, 5, 11, 31, 161, 686]
assert [int((sx <= 10**k).sum()) for k in range(1, 9)] == [1, 2, 2, 2, 5, 5, 18, 82]

# 4020: bi-twin chains below 10^5
bt = [n for n in range(2, 10**5) if isprime(n - 1) and isprime(n + 1) and isprime(2 * n - 1) and isprime(2 * n + 1)]
assert bt[:10] == [6, 30, 660, 810, 2130, 2550, 3330, 3390, 5850, 6270] and len(bt) == 29 and all(n % 30 == 0 for n in bt[1:])

# 4021: Nelson's magic square of nine consecutive primes
c_ = 1480028171; sq = [[c_ + 12, c_ + 18, c_ - 30], [c_ - 42, c_, c_ + 42], [c_ + 30, c_ - 18, c_ - 12]]
assert sorted(x for r in sq for x in r) == [c_ + d for d in (-42, -30, -18, -12, 0, 12, 18, 30, 42)]
assert all(isprime(x) for r in sq for x in r) and nextprime(c_ - 43) == c_ - 42 and nextprime(c_ + 42) > c_ + 42 and not any(isprime(x) for x in range(c_ - 41, c_ + 42) if x - c_ not in (-42, -30, -18, -12, 0, 12, 18, 30, 42))
assert len({sum(r) for r in sq} | {sum(col) for col in zip(*sq)} | {sq[0][0] + sq[1][1] + sq[2][2], sq[0][2] + sq[1][1] + sq[2][0]}) == 1 and sum(sq[0]) == 4440084513

# 4021: the square from the progression 199 + 210i has magic constant 3117
assert all(isprime(199 + 210 * i) for i in range(9)) and 3 * (199 + 210 * 4) == 3117

# 2037: r(N) reaches its trivial upper bound only for N <= 210 (N <= 2e5), via an FFT count
ev = np.zeros(200_001); ev[P[(P > 2) & (P <= 200_000)]] = 1
rN = np.rint(np.fft.irfft(np.fft.rfft(ev, 2**19) ** 2, 2**19)[:200_001]).astype(np.int64)
NN = np.arange(4, 200_001, 2); r_part = (rN[NN] + np.where(at(NN // 2) & (NN // 2 > 2), 1, 0)) // 2
r_part[0] = 1                                      # 4 = 2 + 2
ub = np.searchsorted(P, NN - 2, 'right') - np.searchsorted(P, NN // 2 - 1, 'right')
assert NN[r_part == ub].tolist() == [4, 6, 8, 10, 12, 14, 16, 18, 24, 30, 36, 42, 48, 60, 90, 210]

# 1214: Kempner's function, S(n) = n iff n is prime (n > 4), its exceptions, and Tutescu (n <= 10^5)
def kempner(n):
    best = 1
    for q, e in factorint(n).items():
        m = q
        while sum(m // q**i for i in range(1, 40)) < e: m += q
        best = max(best, m)
    return best
K_ = [kempner(n) for n in range(1, 10**5 + 1)]
assert K_[:11] == [1, 2, 3, 4, 5, 3, 7, 4, 6, 5, 11] and all((K_[n - 1] == n) == isprime(n) for n in range(5, 10**5 + 1))
assert [n for n in range(2, 51) if K_[n - 1] != max(factorint(n))] == [4, 8, 9, 12, 16, 18, 24, 25, 27, 32, 36, 45, 48, 49, 50]
assert all(K_[i] != K_[i + 1] for i in range(len(K_) - 1))

# 1519: C(2n, n) is squarefree only for n = 1, 2, 4 (n <= 3000), and v_2(C(2n, n)) = number of 1s in n
def carries_ok(n):                                 # squarefree iff no prime p carries twice adding n + n
    for p in primerange(2, 2 * n + 1):
        if p * p > 2 * n: break
        c = car = 0; m = n
        while m: d = m % p; car = 1 if 2 * d + car >= p else 0; c += car; m //= p
        if c > 1: return False
    return True
assert [n for n in range(1, 3001) if carries_ok(n)] == [1, 2, 4] and all(bin(n).count('1') == (lambda v: (v & -v).bit_length() - 1)(math.comb(2 * n, n)) for n in range(1, 2001))

# 1519: Sylvester-Schur, C(n, k) has a prime factor > k for n >= 2k (n <= 400)
assert all(max(factorint(math.comb(n, k))) > k for n in range(2, 401) for k in range(1, n // 2 + 1))

# 3022: Jacobsthal j(p#) for p = 2..13
def jac(m):
    run = best = 0
    for x in range(1, 2 * m + 1):
        run = run + 1 if math.gcd(x, m) > 1 else 0; best = max(best, run)
    return best + 1
assert [jac(math.prod(primerange(2, q + 1))) for q in (2, 3, 5, 7, 11, 13)] == [2, 4, 6, 10, 14, 22]

# 3022: 16 is an Erdos-Woods number with 2184 and 2200
assert all(math.gcd(x, 2184) > 1 or math.gcd(x, 2200) > 1 for x in range(2185, 2200))

# 6222: multiplicative persistence, the listed values and the least numbers and primes
def pers(n):
    c = 0
    while n >= 10: n = math.prod(map(int, str(n))); c += 1
    return c
assert [pers(x) for x in [0, 10, 25, 39, 77, 679, 6788, 68889, 2677889, 26888999, 3778888999, 277777788888899]] == list(range(12))
assert [next(n for n in range(10**5) if pers(n) == k) for k in range(7)] == [0, 10, 25, 39, 77, 679, 6788]
assert [next(int(p) for p in P[:20000] if pers(int(p)) == k) for k in range(1, 8)] == [11, 29, 47, 277, 769, 8867, 186889]

# --- fourth research pass: computation records (October 2026) ----------------------------------------

# 6108: RSA-260 (862 bits, 2026) and RSA-896 (896 bits, 2026) split into two PRP factors of 130 and 135 digits
for p_, q_, bits, digits in ((4397328654844826923795068102505872571721883526553349659561256924505973939597593482272505698004801207988043088656411102133523080581, 5028695206842569864686141618253083416610081090075366674776775706538324961364412200138116378509733307971876652984898985905923678379, 862, 260), (636606729769440499166579950236036751749912014371509557713570027508971809534551913252252094954941974952859310861988904737359709200557919, 647218161102195448058768698177623951380616936266986989243011933572862870905830904361851542450154852431416136790787107595965374752513489, 896, 270)):
    N_ = p_ * q_
    assert N_.bit_length() == bits and len(str(N_)) == digits and gmpy2.is_prime(p_) and gmpy2.is_prime(q_)

# 6108: the ECM record (83 digits, 7^337 + 1), the p-1 record (66 digits, 960^119 - 1), the p+1 record (60 digits, L_2366)
ecm = 16559819925107279963180573885975861071762981898238616724384425798932514688349020287
pm1 = 672038771836751227845696565342450315062141551559473564642434674541
pp1 = 725516237739635905037132916171116034279215026146021770250523
assert len(str(ecm)) == 83 and pow(7, 337, ecm) == ecm - 1 and gmpy2.is_prime(ecm)
assert len(str(pm1)) == 66 and pow(960, 119, pm1) == 1 and max(factorint(pm1 - 1)) == 2110402817
assert len(str(pp1)) == 60 and int(lucas(2366)) % pp1 == 0

# 3023: the largest known maximal gap, 1854 after 101412319996363309069 (both ends PRP)
gq = 101412319996363309069
assert int(gmpy2.next_prime(gq)) - gq == 1854 and round(1854 / math.log(gq), 2) == 40.25

# --- facts mined from the puzzle archive (October 2026) -------------------------------------------------

# 10114: consecutive p < q < r with p | qr + 1 below 2e8
p3, q3, r3 = P[:-2].astype(np.int64), P[1:-1].astype(np.int64), P[2:].astype(np.int64)
assert p3[((q3 % p3) * (r3 % p3) + 1) % p3 == 0].tolist() == [2, 3, 61]

# 10115: 4p^2 + 1 = 197 (mod 480) for p > 5, and x -> 4x^2 + 1 hits 0 mod 13 within 5 steps
assert {(4 * p * p + 1) % 480 for p in ps if p > 5 and isprime(4 * p * p + 1)} == {197}
def steps13(x):
    k = 0
    while x % 13: x = (4 * x * x + 1) % 13; k += 1
    return k
assert max(steps13(x) for x in range(13)) == 5
ch = [197072563]
while gmpy2.is_prime(4 * ch[-1] ** 2 + 1): ch.append(4 * ch[-1] ** 2 + 1)
assert len(ch) == 5

# 10116: q = 2p + 1 divides R_p exactly when p = 1, 13, 19 (mod 20)
assert all((pow(10, p, 2 * p + 1) == 1) == (p % 20 in (1, 13, 19)) for p in primerange(7, 10**5) if isprime(2 * p + 1))

# 10117: prime factors of gcd(2^p + 1, 3^p + 1) are 1 (mod 6p) and 1, 3 (mod 8); the first composite gcd at p = 2243399
for p in primerange(5, 3000):
    for q in factorint(math.gcd(2**p + 1, 3**p + 1)):
        assert q % (6 * p) == 1 and q % 8 in (1, 3)
p = 2243399; gg = math.gcd(pow(2, p, 53841577 * 174985123) + 1, pow(3, p, 53841577 * 174985123) + 1)
assert gg == 53841577 * 174985123 == (24 * p + 1) * (78 * p + 1) and isprime(24 * p + 1) and isprime(78 * p + 1)

# 10118: Rotondo's Sophie Germain test for 1 < n <= 2000, and the least bases with a^(p-1) = 1 (mod p^2)
fct = 1
for n in range(2, 2001):
    fct *= n - 1
    assert ((fct * fct + 6 * n) % (n * (2 * n + 1)) == 1) == (n % 2 == 1 and isprime(n) and isprime(2 * n + 1))
assert [next(a for a in range(2, p * p) if pow(a, p - 1, p * p) == 1) for p in primerange(3, 32)] == [8, 7, 18, 3, 19, 38, 28, 28, 14, 115]

# 10210: prime-free gaps between consecutive perfect powers below 2e8, and the pair 503284356-503284375
pw2 = sorted({b ** e for e in range(2, 28) for b in range(2, int(round(2e8 ** (1 / e))) + 2) if b ** e <= 2 * 10**8})
assert [(x, y) for x, y in zip(pw2, pw2[1:]) if np.searchsorted(P, y) == np.searchsorted(P, x, 'right')] == [(8, 9), (25, 27), (32, 36), (121, 125), (2187, 2197), (3125, 3136), (32761, 32768), (79507, 79524), (97336, 97344)]
assert 503284356 == 22434**2 and 503284375 == 55**5 and nextprime(503284356) > 503284375

# 10211: nine primes in (32k, 32k + 32) for k = 23175165755516618
k32 = 23175165755516618
assert sum(isprime(x) for x in range(32 * k32, 32 * k32 + 32)) == 9

# 10306: consecutive primes p with 4p^2 + 1 prime: the least run of 5 starts at 1627; a run of 6 at 1648002773
fp = np.array([gmpy2.is_prime(4 * int(p) ** 2 + 1) for p in P[:20000]])
run5 = next(i for i in range(len(fp) - 5) if fp[i:i + 5].all())
assert int(P[run5]) == 1627
r6 = [1648002773]
while len(r6) < 7: r6.append(nextprime(r6[-1]))
assert [gmpy2.is_prime(4 * x * x + 1) for x in r6] == [True] * 6 + [False] and not gmpy2.is_prime(4 * int(gmpy2.prev_prime(1648002773)) ** 2 + 1)

# 10307: 2^n and F_n as sums of two consecutive primes
def two_consec(N):                                  # N = a + b with a < b consecutive primes
    a_ = gmpy2.prev_prime(N // 2 + 1); return a_ + gmpy2.next_prime(a_) == N
assert [n for n in range(3, 401) if two_consec(2**n)] == [3, 7, 13, 77, 182]
assert [n for n in range(5, 301) if two_consec(int(fibonacci(n)))] == [5, 6, 12, 60]

# 10308: the run of 15 integers with Omega = 4 from 97524222465 is exact
om4 = lambda n: sum(factorint(n).values())
assert all(om4(97524222465 + i) == 4 for i in range(15)) and om4(97524222464) != 4 and om4(97524222480) != 4

# 10406: twin chains (p, p + 2) -> (2p + 1, 2p + 3): least starts for k = 1..4, and the 5- and 6-chains
def tchain(p):
    k = 0
    while isprime(p) and isprime(p + 2): k += 1; p = 2 * p + 1
    return k
assert [next(int(p) for p in P[:30000] if tchain(int(p)) >= k) for k in (1, 2, 3, 4)] == [3, 5, 211049, 253679]
assert tchain(41887255409) == 5 and tchain(73768891456259) == 6

# 10607: p, 2p + 1, ..., np + (n - 1) all prime: least p for n = 6, 7, 8; the chains for n = 9, 10
def kchain(p):
    k = 1
    while gmpy2.is_prime((p + 1) * (k + 1) - 1): k += 1
    return k
assert [next(int(p) for p in P[:300000] if kchain(int(p)) >= n) for n in (6, 7, 8)] == [154769, 2894219, 2894219]
assert kchain(407874179) >= 9 and kchain(214580145779) >= 10

# 10607: least m with 2^(2^k) + m prime for k = 0..n (n = 5, 6), and the values for n = 7, 8
fm_ = lambda m, n: all(gmpy2.is_prime(2**(2**k) + m) for k in range(n + 1))
assert next(m for m in range(1, 10**5) if fm_(m, 5)) == 15 and next(m for m in range(1, 10**5) if fm_(m, 6)) == 66747
assert fm_(475425, 7) and fm_(12124167, 8)

# 10608: n! + 239 prime list, 648k^6 + 1, floor(n^2/3)
f, l239 = 1, []
for n in range(1, 239):
    f *= n
    if gmpy2.is_prime(f + 239): l239.append(n)
assert l239 == [2, 4, 5, 7, 8, 9, 10, 11, 18, 22, 38, 45, 48, 76, 90, 96, 99, 132, 136, 147, 225]
assert next(648 * k**6 + 1 for k in range(1, 100) if isprime(648 * k**6 + 1)) == 472393
assert not any(isprime(n * n // 3) for n in range(5, 20000))

# 10706: weakly primes below 10^6, and the period-666 prime
def weak(p):
    s_ = str(p)
    return not any(isp[int(s_[:i] + d + s_[i + 1:])] for i in range(len(s_)) for d in '0123456789' if d != s_[i])
assert [p for p in ps + [int(x) for x in P[(P >= 10**5) & (P < 10**6)]] if weak(p)][:5] == [294001, 505447, 584141, 604171, 971767]
assert n_order(10, 902659997773) == 666 and isprime(902659997773)

# 10707: Chernick numbers with one prime factor (m < 10^5), and a Carmichael number 2n + 1
def korselt(N):
    fs = factorint(N); return len(fs) > 1 and all(e == 1 for e in fs.values()) and all((N - 1) % (q - 1) == 0 for q in fs)
one = [m for m in range(1, 10**5) if pow(2, (6*m+1)*(12*m+1)*(18*m+1) - 1, (6*m+1)*(12*m+1)*(18*m+1)) == 1
       and sum(isprime(x) for x in (6*m+1, 12*m+1, 18*m+1)) == 1 and korselt((6*m+1)*(12*m+1)*(18*m+1))]
assert one == [22, 105, 225, 609, 12259, 15125, 37961, 85410]
assert 82163 * 1068107 * 3532967 == 310049210890163447 and korselt(310049210890163447) and not korselt(310049210890163447 // 2)

# 10708: covers for the Cohen-Selfridge prime and the prime Brier number; Cullen identities
kcs, brier = 47867742232066880047611079, 10439679896374780276373
cov48, cov360 = (3, 5, 13, 17, 97, 241, 257), (3, 7, 11, 19, 31, 37, 41, 61, 73, 109, 151, 331)
assert gmpy2.is_prime(kcs) and gmpy2.is_prime(brier)
assert all(any((kcs + pow(2, n, q)) % q == 0 for q in cov48) and any((kcs - pow(2, n, q)) % q == 0 for q in cov360) for n in range(720))
assert all(any((brier * pow(2, n, q) + 1) % q == 0 for q in cov48) and any((brier * pow(2, n, q) - 1) % q == 0 for q in cov360) for n in range(720))
assert all(n * 2**(2 * n + 1) + 1 == 2 * n * 2**(2 * n) + 1 for n in range(1, 50))
assert all((p * pow(2, p, 3) + 1) % 3 == 0 for p in ps if p % 6 == 1) and all((p * pow(2, p, 7) + 1) % 7 == 0 for p in ps if p % 42 == 5)

# 10709: emirp counts by digits, and 14771 for the Wilson quotient
emp = [int((at(rev(P[(P >= 10**(d - 1)) & (P < 10**d)])) & (rev(P[(P >= 10**(d - 1)) & (P < 10**d)]) != P[(P >= 10**(d - 1)) & (P < 10**d)])).sum()) for d in range(1, 8)]
assert emp == [0, 8, 28, 204, 1406, 9538, 70474]
p = 14771; wq = (math.factorial(p - 1) + 1) // p
assert (math.factorial(p - 1) + 1) % (p * p) != 0 and pow(wq, p - 1, p * p) == 1

# 10805: de Polignac counts below 10^7 and Erdos's class
odd = np.arange(3, 10**7, 2); hitp = np.zeros(len(odd), bool)
for k in range(24): v_ = odd - 2**k; hitp |= at(np.maximum(v_, 0)) & (v_ > 1)
pol_ = odd[~hitp]
assert len(pol_) == 421862 and int(isp[pol_].sum()) == 102612 and int((~isp[pol_]).sum()) == 319250
assert all(any((18814027 - pow(2, k, q)) % q == 0 for q in (3, 5, 7, 13, 17, 241)) for k in range(24)) and 346729110 == 62 * 5592405

print("all FACTS.md checks passed")
