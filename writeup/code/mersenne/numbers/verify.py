"""Independent recheck (plain Python integers and Fractions, no numpy, no grammar.py).

1. The real sequence is rebuilt from scratch: sympy.isprime on each term (< 2^64, so a proof) and the
   Lucas-Lehmer test for the terms below 5000 (gmpy2-free, plain ints).
2. Sub-grammar G0 re-implemented: x+c (|c|<=30), a*x+d (a in 2..8,10,16; |d|<=15), floor/round/ceil(x*p/q)
   for the rationals of grammar.py, next y > x with y = r mod m, next binary / decimal palindrome,
   binary and decimal reversal, next number with the same popcount, rotl1/rotr1.
   For every one of the 501 sequences the best exact-transition count over G0 is recomputed and compared with
   the numpy grammar on the same names (the numpy run is done here too, so the two must agree row by row).
3. The named top rules (RULES below) are recomputed on the real sequence: hit positions and the run from a_1.
"""
import sys, math, pickle
from fractions import Fraction
import numpy as np
import sympy
sys.path.insert(0, '.')

D = np.load('seqs.npz')
rows = [list(map(int, D['real']))] + [list(map(int, r)) for r in D['nullA']] + [list(map(int, r)) for r in D['nullB']]
real = rows[0]

# 1. the sequence
def ll(p):
    if p == 2: return True
    M = (1 << p) - 1; s = 4
    for _ in range(p - 2): s = (s * s - 2) % M
    return s == 0
assert all(sympy.isprime(v) for v in real)
ll_small = [p for p in sympy.primerange(2, 5000) if ll(p)]
print('1. all 52 terms prime (sympy, proof < 2^64):', True, '| Lucas-Lehmer below 5000 gives', ll_small == [v for v in real if v < 5000])

# 2. sub-grammar
def bl(x): return x.bit_length()
def revb(x): return int(bin(x)[2:][::-1], 2)
def revd(x): return int(str(x)[::-1])
def isbpal(x): s = bin(x)[2:]; return s == s[::-1]
def isdpal(x): s = str(x); return s == s[::-1]
def nextp(x, pred):
    base = 2 if pred is isbpal else 10
    y = x + 1
    d = []
    t = y
    while t: d.append(t % base); t //= base
    d = d[::-1]; L = len(d); h = (L + 1) // 2
    def mirror(first):
        full = first + first[:L // 2][::-1]
        v = 0
        for c in full: v = v * base + c
        return v
    c = mirror(d[:h])
    if c >= y: return c
    v = 0
    for q in d[:h]: v = v * base + q
    v += 1
    if v == base ** h:  # carry into a longer length: base^L + 1
        return base ** L + 1
    first = []
    for _ in range(h): first.append(v % base); v //= base
    return mirror(first[::-1])
def gos(x):
    c = x & -x; r = x + c; return (((r ^ x) >> 2) // c) | r
def rotl(x): L = bl(x); return ((x << 1) & ((1 << L) - 1)) | (x >> (L - 1))
def rotr(x): L = bl(x); return (x >> 1) | ((x & 1) << (L - 1))

G0 = {}
for c in range(-30, 31):
    if c: G0[f'x{c:+d}'] = (lambda c: lambda x: x + c)(c)
for a in (2, 3, 4, 5, 6, 7, 8, 10, 16):
    for d in range(-15, 16):
        G0[f'{a}x{d:+d}'] = (lambda a, d: lambda x: a * x + d)(a, d)
RATS = {}
for q in range(1, 7):
    for p in range(q + 1, 4 * q + 1):
        fr = Fraction(p, q)
        if fr.denominator == q and fr != int(fr): RATS[f'{p}/{q}'] = fr
for nm, fr in list(RATS.items()):
    RATS['1/' + nm] = 1 / fr
for nm, fr in RATS.items():
    G0[f'floor(x*{nm})'] = (lambda fr: lambda x: math.floor(x * fr))(fr)
    G0[f'ceil(x*{nm})'] = (lambda fr: lambda x: math.ceil(x * fr))(fr)
    # numpy rint rounds halves to even; Fraction round() does the same
    G0[f'round(x*{nm})'] = (lambda fr: lambda x: round(x * fr))(fr)
for m in list(range(2, 13)) + [16, 18, 20, 24, 30, 60]:
    for r in range(m):
        G0[f'next≡{r}(mod {m})'] = (lambda m, r: lambda x: x + 1 + ((r - x - 1) % m))(m, r)
G0['revbin'] = revb; G0['revdec'] = revd; G0['nextsamepop'] = gos; G0['rotl1'] = rotl; G0['rotr1'] = rotr
G0['next binpal'] = lambda x: nextp(x, isbpal) if x < (1 << 40) else -1
G0['next decpal'] = lambda x: nextp(x, isdpal) if x < 10 ** 12 else -1


def best_exact(fns, seqs):
    out = []
    for s in seqs:
        b = 0
        for f in fns:
            h = sum(1 for i in range(len(s) - 1) if f(s[i]) == s[i + 1])
            b = max(b, h)
        out.append(b)
    return out

# the palindrome successors are slow in plain Python for big gaps; they are checked on every row anyway
plain = best_exact(list(G0.values()), rows)
# numpy grammar on the same names
from grammar import build, apply
P = {nm: fn for nm, fn, _ in build()}
missing = [nm for nm in G0 if nm not in P]
M = np.array(rows, np.int64); X = M[:, :-1]; Y = M[:, 1:]; N = np.broadcast_to(np.arange(1, 52), X.shape).astype(np.int64)
nb = np.zeros(len(rows), np.int64)
for nm in G0:
    if nm in P:
        nb = np.maximum(nb, (apply(P[nm], X, N) == Y).sum(1))
diff = [i for i in range(len(rows)) if plain[i] != nb[i]]
print(f'2. sub-grammar G0: {len(G0)} maps (names missing in grammar.py: {missing}); best exact per row, plain vs numpy: '
      f'{len(rows) - len(diff)} of {len(rows)} rows agree; differing rows: {diff[:10]}')
print(f'   real {plain[0]}; nullA median {sorted(plain[1:251])[125]}; nullB median {sorted(plain[251:])[125]}; '
      f'p(nullB >= real) = {(1 + sum(v >= plain[0] for v in plain[251:])) / 251:.3f}')

# 3. named rules on the real sequence
def nbp(x): return nextp(x, isbpal)
def pm1(x):
    return min(v for k in range(1, 62) for v in ((1 << k) - 1, (1 << k) + 1) if v > x)
RULES = {
    'revbin(x + (x>>1) | 1)': lambda x: revb((x + (x >> 1)) | 1),
    'next binpal(floor(7x/6))': lambda x: nbp(math.floor(Fraction(7, 6) * x)),
    'next 2^k+-1 after floor(7x/6)': lambda x: pm1(math.floor(Fraction(7, 6) * x)),
    'x-2 o next2^k+-1 o next2^k+-1': lambda x: pm1(pm1(x)) - 2,
    'rotl1(x + revbin(x))': lambda x: rotl(x + revb(x)),
    'revbin(x) + 30': lambda x: revb(x) + 30,
    'ceil((3x+n)/phi^2) [n = index]': None,
    'x+2': lambda x: x + 2,
    '2x-1': lambda x: 2 * x - 1,
}
extra = {}
try:
    extra = pickle.load(open('verify_rules.pkl', 'rb'))  # optional: (name, python-lambda-source) pairs added later
except Exception:
    pass
for nm, src in extra.items():
    RULES[nm] = eval(src)
for nm, f in RULES.items():
    if f is None: continue
    hits = [i for i in range(51) if f(real[i]) == real[i + 1]]
    run0 = 0
    while run0 < 51 and f(real[run0]) == real[run0 + 1]: run0 += 1
    print(f'3. {nm:40s} exact {len(hits):2d}, run from a_1 {run0}, hits: ' + ' '.join(f'{real[i]}->{real[i + 1]}' for i in hits))

# 4. floor(c*r^n) with the float parameters found by formulas.py, recomputed exactly with Fractions
c, r = (Fraction(float(t)) for t in open('floor_cr_params.txt').read().split())
hits = [n for n in range(1, 53) if math.floor(c * r ** n) == real[n - 1]]
print(f'4. floor(c*r^n), c = {float(c)!r}, r = {float(r)!r} (exact rational arithmetic): hits at n = {hits}')
# 5. the 3x+n rule with ceil(x / phi^2): phi^2 irrational, use a 60-digit decimal
from decimal import Decimal, getcontext
getcontext().prec = 60
phi2 = ((1 + Decimal(5).sqrt()) / 2) ** 2
h = [n for n in range(1, 52) if int(((3 * real[n - 1] + n) / phi2).to_integral_value(rounding='ROUND_CEILING')) == real[n]]
print('5. ceil((3 a_n + n)/phi^2) = a_(n+1) at n =', h)
