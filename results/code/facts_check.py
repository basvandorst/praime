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
# --- facts added in the October 2026 research pass (sections 1B-1G, 2, 4, 5, 6) ------------------
from sympy import cyclotomic_poly, primefactors, totient
isp = np.zeros(200_000_001, dtype=bool); isp[P] = True
def at(a): return np.where(a <= 200_000_000, isp[np.minimum(a, 200_000_000)], False)
def upto(Q, k): return Q[Q <= 10**k]
tup = {name: [int(f(upto(P, k)).sum()) for k in range(1, 9)] for name, f in {
    'cousin': lambda Q: at(Q + 4), 'sexy': lambda Q: at(Q + 6), 'sg': lambda Q: at(2 * Q + 1),
    'trip': lambda Q: (at(Q + 2) & at(Q + 6)) | (at(Q + 4) & at(Q + 6)), 'quad': lambda Q: at(Q + 2) & at(Q + 6) & at(Q + 8)}.items()}
assert tup['cousin'] == [2, 9, 41, 203, 1216, 8144, 58622, 440258] and tup['sexy'] == [2, 16, 74, 411, 2447, 16386, 117207, 879908]
assert tup['sg'] == [3, 10, 37, 190, 1171, 7746, 56032, 423140] and tup['trip'] == [2, 9, 30, 112, 507, 2837, 17220, 111156]
assert tup['quad'] == [1, 2, 5, 12, 38, 166, 899, 4768]
def cc_least(d):                                                  # Cunningham chains p, 2p + d, ...: least start with >= k primes
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
for p, d in ((75949683268775452969189499, 1), (3352894966466100254851081, -1)):   # length-20 chains (Armitage 2026)
    L, x = 0, p
    while gmpy2.is_prime(x): L += 1; x = 2 * x + d
    assert L == 20 and not gmpy2.is_prime((p - d) // 2)
for m, cl, x0 in ((8, (1, 3, 5, 7), 23), (12, (1, 5, 7, 11), 463)):      # races: the class of 1 is strictly last
    C = np.stack([np.cumsum(P % m == c) for c in cl]); i0 = int(np.searchsorted(P, x0))
    assert np.all(C[0, i0:] < C[1:, i0:].min(0)) and not np.all(C[0, i0 - 1:] < C[1:, i0 - 1:].min(0))
Q = P[P > 5]; C = {c: np.cumsum(Q % 10 == c) for c in (1, 3, 7, 9)}
assert [round(100 * float(np.mean(C[a] > C[b])), 1) for a, b in ((9, 3), (9, 7), (1, 3), (1, 7))] == [1.0, 1.2, 2.5, 2.8]
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
c4 = np.cumsum(P % 4 == 1) - np.cumsum(P % 4 == 3)
assert [int(x) for x in P[c4 == 0][:6]] == [2, 5, 17, 41, 461, 26833]
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
def g0(p):
    qs = primefactors(p - 1); a = 2
    while any(pow(a, (p - 1) // q, p) == 1 for q in qs): a += 1
    return a
assert records(g0, 10**6) == [(2, 3), (3, 7), (5, 23), (6, 41), (7, 71), (19, 191), (21, 409), (23, 2161), (31, 5881), (37, 36721), (38, 55441), (44, 71761), (69, 110881), (73, 760321)]
P7 = [int(p) for p in P[P < 10**7]]
for b, known in ((3, [11, 1006003]), (5, [2, 20771, 40487, 53471161]), (7, [5, 491531]), (10, [3, 487, 56598313]), (11, [71])):
    assert [p for p in P7 if b % p and gmpy2.powmod(b, p - 1, p * p) == 1] == [x for x in known if x < 10**7]
    assert all(pow(b, x - 1, x * x) == 1 and isprime(x) for x in known if x > 2)
assert [n for n in range(2, 2000) if sum(pow(i, n - 1, n) for i in range(1, n)) % n == n - 1] == list(primerange(2, 2000))   # Giuga
exc = set()                                                       # Zsigmondy for a^n - 1, a <= 30, n <= 60
for a in range(2, 31):
    for n in range(2, 61):
        v = int(cyclotomic_poly(n, a))
        for q in primefactors(n):
            while v % q == 0: v //= q
        if v == 1: exc.add((a, n))
assert exc == {(2, 6), (3, 2), (7, 2), (15, 2)}
pw = sorted({b ** e for e in range(2, 41) for b in range(2, int(round(10 ** (12 / e))) + 2) if b ** e < 10**12})
assert [x for x, y in zip(pw, pw[1:]) if y - x == 1] == [8]       # Catalan below 10^12
sy = [2]
for _ in range(7): sy.append(sy[-1] ** 2 - sy[-1] + 1)
assert sy[:6] == [2, 3, 7, 43, 1807, 3263443] and all(math.gcd(a, b) == 1 for i, a in enumerate(sy) for b in sy[:i])
assert all(Fraction(1, 1) - sum(Fraction(1, x) for x in sy[:n]) == Fraction(1, sy[n] - 1) for n in range(1, 8))
assert all(q % 6 != 5 for x in sy[:7] for q in factorint(x))
cnt = [0, 0, 0]                                                   # Chebotarev for x^3 - 2: 3, 1, 0 roots
for p in primerange(5, 10**6):
    cnt[1 if p % 3 == 2 else (0 if pow(2, (p - 1) // 3, p) == 1 else 2)] += 1
assert [round(c / sum(cnt), 4) for c in cnt] == [0.1660, 0.5002, 0.3338]
lin = 0                                                           # least prime in progressions vs (phi(q) ln q)^2
for q in range(4, 301):
    r = P[:200000] % q; u, idx = np.unique(r, return_index=True); ok = np.gcd(u, q) == 1
    assert ok.sum() == totient(q)
    lin = max(lin, int(P[idx[ok]].max()) / (int(totient(q)) * math.log(q)) ** 2)
assert 0.650 < lin < 0.651
assert all(2 ** d - 1 >= min(int(p) for p in P[:50000] if p % d == 1) for d in range(3, 200))
irr = {}                                                          # irregular primes via p B_k = sum j^k (mod p^2)
for p in primerange(5, 2000):
    m = p * p; j = np.arange(1, p, dtype=np.int64); j2 = j * j % m; pwr = j2.copy(); ks = []
    for k in range(2, p - 2, 2):
        if int(pwr.sum() % m) == 0: ks.append(k)
        pwr = pwr * j2 % m
    if ks: irr[p] = ks
assert list(irr)[:12] == [37, 59, 67, 101, 103, 131, 149, 157, 233, 257, 263, 271] and sum(p < 1000 for p in irr) == 64
assert len(irr) == 121 and len(list(primerange(5, 2000))) == 301 and irr[157] == [62, 110] and irr[491] == [292, 336, 338]
assert min(p for p in irr if len(irr[p]) == 3) == 491 and min(p for p in irr if len(irr[p]) == 2) == 157
p = 16843; m = p**4; c = 1
for i in range(1, p): c = c * (p + i) % m * pow(i, -1, m) % m     # C(2p-1, p-1) = prod (p+i)/i
assert c == 1
# section 2: Mertens, Liouville, Robin, Lagarias, omega, Goldbach
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
Y = 10**6; sig = np.zeros(Y + 1, np.int64)
for d in range(1, Y + 1): sig[d::d] += d
nn = np.arange(3, Y + 1); viol = nn[sig[3:] >= np.exp(0.5772156649015329) * nn * np.log(np.log(nn))]
assert viol.tolist() == [3, 4, 5, 6, 8, 9, 10, 12, 16, 18, 20, 24, 30, 36, 48, 60, 72, 84, 120, 180, 240, 360, 720, 840, 2520, 5040]
assert round(19344 / (math.exp(0.5772156649015329) * 5040 * math.log(math.log(5040))), 4) == 1.0056 and sig[5040] == 19344
H = np.cumsum(1 / np.arange(1, 10**5 + 1)); n5 = np.arange(2, 10**5 + 1)
assert np.all(sig[2:10**5 + 1] < H[1:] + np.exp(H[1:]) * np.log(H[1:]))           # Lagarias
phi = np.arange(Y + 1, dtype=np.int64)
for p in P[P <= Y].tolist(): phi[::p] -= phi[::p] // p
comp = np.nonzero(~isp[:Y + 1])[0]; comp = comp[comp > 1]
assert not np.any((comp - 1) % phi[comp] == 0)                                     # Lehmer's totient problem
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
gn = 3325581707333960528
assert isprime(gn - 9781) and not any(isprime(q) and isprime(gn - q) for q in range(3, 9781, 2))
# section 5 and 6
c2, out, nxt = 0, [], 10
for n in range(1, 10**6 + 1):
    c2 += gmpy2.is_prime(n * n + 1)
    if n == nxt: out.append(c2); nxt *= 10
assert out == [5, 19, 112, 841, 6656, 54110]
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
assert min(n for n in sp if gmpy2.is_strong_prp(n, 3)) == 1373653 and not any(gmpy2.is_strong_prp(n, 3) and gmpy2.is_strong_prp(n, 5) for n in sp)
car = [int(n) for n in psp if all(e == 1 for e in factorint(int(n)).values()) and all((n - 1) % (q - 1) == 0 for q in factorint(int(n)))]
assert [sum(x <= 10**k for x in car) for k in range(3, 8)] == [1, 7, 16, 43, 105]
mr = [2047, 1373653, 25326001, 3215031751, 2152302898747, 3474749660383, 341550071728321, 341550071728321,
      3825123056546413051, 3825123056546413051, 3825123056546413051, 318665857834031151167461, 3317044064679887385961981]
for n, x in enumerate(mr, 1):
    assert not isprime(x) and all(gmpy2.is_strong_prp(x, b) for b in list(primerange(2, 50))[:n])
    assert n == 13 or not all(gmpy2.is_strong_prp(x, b) for b in list(primerange(2, 50))[:n + 1]) or mr[n] == x
for k in range(1, 100):                                                             # Chernick
    a, b, c = 6 * k + 1, 12 * k + 1, 18 * k + 1
    if isprime(a) and isprime(b) and isprime(c): assert all((a * b * c - 1) % (q - 1) == 0 for q in (a, b, c))
assert all((pow(2, p, 2 * p + 1) == 1) == isprime(2 * p + 1) for p in primerange(5, 10**5) if p % 4 == 3)   # Euler-Lagrange
for k, cover in ((509203, (3, 5, 7, 13, 17, 241)), (271129, (3, 5, 7, 13, 17, 241))):
    sgn = -1 if k == 509203 else 1
    assert all(any((k * pow(2, n, q) + sgn) % q == 0 for q in cover) for n in range(1, 200))
f, fp, fm = 1, [], []
for n in range(1, 501):
    f *= n
    if gmpy2.is_prime(f + 1): fp.append(n)
    if gmpy2.is_prime(f - 1): fm.append(n)
assert fp == [1, 2, 3, 11, 27, 37, 41, 73, 77, 116, 154, 320, 340, 399, 427]
assert fm == [3, 4, 6, 7, 12, 14, 30, 32, 33, 38, 94, 166, 324, 379, 469]
pr, pp, pm = 1, [], []
for p in primerange(2, 2000):
    pr *= p
    if gmpy2.is_prime(pr + 1): pp.append(p)
    if gmpy2.is_prime(pr - 1): pm.append(p)
assert pp == [2, 3, 5, 7, 11, 31, 379, 1019, 1021] and pm == [3, 5, 11, 13, 41, 89, 317, 337, 991, 1873]
pr, fo = 1, []
for p in primerange(2, 600):
    pr *= p; fo.append(int(gmpy2.next_prime(pr + 1)) - pr)
assert fo[:10] == [3, 5, 7, 13, 23, 17, 19, 23, 37, 61] and all(isprime(x) for x in fo[:100])
f = 1; assert [n for n in range(1, 1001) if gmpy2.is_square((f := f * n) + 1)] == [4, 5, 7]   # Brocard
em = [2]
while len(em) < 9: em.append(min(factorint(math.prod(em) + 1)))
assert em == [2, 3, 7, 43, 13, 53, 5, 6221671, 38709183810571]
pier = sorted(2**u * 3**v + 1 for u in range(22) for v in range(14) if 2**u * 3**v + 1 < 2 * 10**6 and isprime(2**u * 3**v + 1))
assert len(pier) == 46 and pier[:12] == [2, 3, 5, 7, 13, 17, 19, 37, 73, 97, 109, 163]
assert [p for p in primerange(3, 1000) if gmpy2.is_prime((2**p + 1) // 3)] == [3, 5, 7, 11, 13, 17, 19, 23, 31, 43, 61, 79, 101, 127, 167, 191, 199, 313, 347, 701]
assert [n for n in range(2, 1100) if gmpy2.is_prime((10**n - 1) // 9)] == [2, 19, 23, 317, 1031]
assert [n for n in range(1, 1001) if gmpy2.is_prime(n * 2**n + 1)] == [1, 141] and [n for n in range(1, 1001) if gmpy2.is_prime(n * 2**n - 1)] == [2, 3, 6, 30, 75, 81, 115, 123, 249, 362, 384, 462, 512, 751, 822]
def pal(d):
    h = (d + 1) // 2; out = 0
    for x in range(10**(h - 1), 10**h):
        s_ = str(x); v = int(s_ + s_[::-1][d % 2:])
        out += bool(isp[v]) if v <= 200_000_000 else isprime(v)
    return out
assert [pal(d) for d in range(1, 10)] == [4, 1, 15, 0, 93, 0, 668, 0, 5172]
def rev(a):
    r = np.zeros_like(a); a = a.copy()
    while a.any(): r = r * 10 + a % 10; a //= 10
    return r
assert [int(at(rev(P[(P >= 10**(d - 1)) & (P < 10**d)])).sum()) for d in range(1, 8)] == [4, 9, 43, 204, 1499, 9538, 71142]
assert [k for k in range(101) if gmpy2.is_prime(10**(2 * k + 4) + 666 * 10**(k + 1) + 1)] == [0, 13, 42]
gp = lambda R: sum(1 for a in range(-R, R + 1) for b in range(-R, R + 1) if a * a + b * b <= R * R and
                   (isprime(a * a + b * b) if a and b else (isprime(abs(a + b)) and abs(a + b) % 4 == 3)))
assert gp(10) == 100 and gp(100) == 4928
assert int(gmpy2.next_prime(10**999)) == 10**999 + 7
print("all FACTS.md checks passed")
