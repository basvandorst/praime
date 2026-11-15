"""Independent recheck of the real-sequence statistics with plain integers / Fractions
and different algorithms (no numpy for the arithmetic, no code shared with angle*.py)."""
import json, math, sys
from fractions import Fraction
from itertools import product
import sympy

A = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217,
     4253, 4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243, 110503, 132049,
     216091, 756839, 859433, 1257787, 1398269, 2976221, 3021377, 6972593, 13466917,
     20996011, 24036583, 25964951, 30402457, 32582657, 37156667, 42643801, 43112609,
     57885161, 74207281, 77232917, 82589933, 136279841]
assert all(sympy.isprime(x) for x in A)
n = len(A)
res = {}

def rank_int(M):
    """exact rank by fraction-free (Bareiss-style) elimination over the integers"""
    M = [list(r) for r in M]
    rk = 0
    cols = len(M[0]) if M else 0
    prev = 1
    for c in range(cols):
        piv = next((r for r in range(rk, len(M)) if M[r][c] != 0), None)
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        for r in range(rk + 1, len(M)):
            row = []
            for j in range(cols):
                q, rem = divmod(M[rk][c] * M[r][j] - M[r][c] * M[rk][j], prev)
                assert rem == 0
                row.append(q)
            M[r] = row
        prev = M[rk][c]
        rk += 1
    return rk

# 1. linear complexity via ranks of Hankel systems over Q (sympy), m = 1..52
def lin_cx(m, start=0):
    # linear complexity is nondecreasing in m, so the search may start at L(m-1)
    for L in range(start, m + 1):
        if L == 0:
            if all(x == 0 for x in A[:m]):
                return 0
            continue
        rows = [[A[t - i] for i in range(1, L + 1)] for t in range(L, m)]
        rhs = [A[t] for t in range(L, m)]
        if not rows:
            return L
        Mb = [r + [v] for r, v in zip(rows, rhs)]
        if rank_int(rows) == rank_int(Mb):
            return L
    return m

Ls = []
for m in range(1, n + 1):
    Ls.append(lin_cx(m, Ls[-1] if Ls else 0))
# zero discrepancy at term m+1 <=> L(m+1) == L(m) and the recurrence of length L(m) is unique/extends:
# use the definition directly: some recurrence of length L(m) fitting m terms also fits m+1 terms
res['bm_def'] = max(math.ceil(m / 2) - Ls[m - 1] for m in range(8, n + 1))
res['linear complexity profile m=1..52'] = Ls
# unique shortest recurrence (2L(m) <= m) that also fits term m+1
res['bm_zero(unique, L(m+1)==L(m), m>=8)'] = sum(Ls[m] == Ls[m - 1] and 2 * Ls[m - 1] <= m for m in range(8, n))

# 1b. 2,3,5,7,13: a_{n+2} = 4 a_n - a_{n+1}
res['4*5-7==13'] = 4 * 5 - 7 == 13
res['next 4*7-13'] = 4 * 7 - 13

# 1c. small integer recurrences, plain loops
def runs_small(const):
    best = 0
    lead = 0
    for k in range(1, 5):
        for c in product(range(-10, 11), repeat=k):
            R = [A[t] - sum(c[i] * A[t - 1 - i] for i in range(k)) for t in range(k, n)]
            if const:
                cur = 1
                b = 1
                for t in range(1, len(R)):
                    cur = cur + 1 if R[t] == R[t - 1] else 1
                    b = max(b, cur)
                l = 1
                while l < len(R) and R[l] == R[0]:
                    l += 1
            else:
                cur = b = 0
                for r in R:
                    cur = cur + 1 if r == 0 else 0
                    b = max(b, cur)
                l = 0
                while l < len(R) and R[l] == 0:
                    l += 1
            best = max(best, b)
            lead = max(lead, l)
    return best, lead

res['si_run_h, si_lead_h'] = runs_small(False)
res['si_run_c, si_lead_c'] = runs_small(True)

# 2. relations, plain loops
S = {x: i for i, x in enumerate(A)}
for D in (0, 1, 2, 10):
    tot = 0
    tot1 = 0
    has = set()
    for k in range(n):
        for j in range(k):
            for i in range(j):
                for c1 in range(-8, 9):
                    if c1 == 0:
                        continue
                    for c2 in range(-8, 9):
                        if c2 == 0:
                            continue
                        if abs(A[k] - c1 * A[i] - c2 * A[j]) <= D:
                            tot += 1
                            has.add(k)
            for c in range(1, 9):
                if abs(A[k] - c * A[j]) <= D:
                    tot1 += 1
                    has.add(k)
    res['L2tot_D%d' % D] = tot
    res['L1tot_D%d' % D] = tot1
    res['Lhas_D%d' % D] = sum(1 for k in has if k >= 8)

mul = 0
mulhas = set()
for k in range(n):
    for i in range(k):
        for j in range(i, k):
            for l in range(k):
                if l in (i, j):
                    continue
                q, r = divmod(A[i] * A[j], A[l])
                v = q + (1 if 2 * r >= A[l] else 0)
                if v == A[k]:
                    mul += 1
                    mulhas.add(k)
res['mul_tot'] = mul
res['mul_has'] = sum(1 for k in mulhas if k >= 8)

pw = 0
for q in range(1, 5):
    for p in range(q + 1, 8 * q + 1):
        if math.gcd(p, q) != 1:
            continue
        for i in range(n):
            for k in range(i + 1, n):
                # A[k] == round(A[i]^(p/q))  <=>  (2A[k]-1)^q <= 2^q A[i]^p < (2A[k]+1)^q
                lhs = (2 * A[k] - 1) ** q
                mid = 2 ** q * A[i] ** p
                rhs = (2 * A[k] + 1) ** q
                if lhs <= mid < rhs:
                    pw += 1
res['pow_exact'] = pw

# 3. bits: LZ76 by the phrase definition; min DFA by Myhill-Nerode right languages
def lz76_naive(s):
    c = 0
    i = 0
    N = len(s)
    while i < N:
        L = 1
        # extend while s[i:i+L] occurs in s[0:i+L-1]
        while i + L <= N and s[i:i + L] in s[0:i + L - 1]:
            L += 1
        c += 1
        i += L
    return c

cat = ''.join(bin(x)[2:] for x in A)
res['lz_cat'] = lz76_naive(cat) * math.log2(len(cat)) / len(cat)

def nerode(strings):
    Lset = set(strings)
    prefixes = {w[:i] for w in Lset for i in range(len(w) + 1)}
    classes = {frozenset(s[len(w):] for s in Lset if s.startswith(w)) for w in prefixes}
    return len(classes)

tv = sv = 0
for k in range(8, 21):
    sub = [bin(x)[2:] for x in A if x < 2 ** k]
    sv += nerode(sub)
    tv += len({w[:i] for w in sub for i in range(len(w) + 1)})
res['dfa_var'] = sv / tv

# 5. growth
x = [math.log2(v) for v in A]
g = [x[i + 1] - x[i] for i in range(n - 1)]
mu = sum(g) / len(g)
sd = math.sqrt(sum((v - mu) ** 2 for v in g) / len(g))
res['gap_cv'] = sd / mu
res['gap_small'] = sum(1 for v in g if v < 0.1 * mu)
res['rate terms/bit (51 gaps / span)'] = 51 / (x[-1] - x[0])

json.dump(res, open('check_real.json', 'w'), indent=1, default=str)


for k, v in res.items():
    print(k, ':', v, flush=True)
