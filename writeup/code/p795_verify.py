# Puzzle 795: F(a, m) reverses the digit strings between successive occurrences of the digit a in m.
# Checks (1) every solution in data/p795_solutions.txt: p prime, exactly D distinct digits, every F(a, p) a prime
# different from p (sympy.isprime is a proof below 2^64; PARI isprime for the larger ones), (2) Vantieghem's 18- and
# 21-digit numbers, (3) the lemma that every used digit occurs at least twice and that 10 digits need 21 digits,
# (4) brute force for L = 7, 8 (all numbers with D distinct digits, each at least twice) as an independent check
# of the C search.   Run from the repo root: .venv/bin/python3 writeup/code/p795_verify.py
import itertools
from sympy import isprime
from sympy.utilities.iterables import multiset_permutations
import cypari2
pari = cypari2.Pari()

def F(a, s):
    pos = [i for i, c in enumerate(s) if c == a]
    t = list(s)
    for i, j in zip(pos, pos[1:]):
        t[i + 1:j] = t[i + 1:j][::-1]
    return "".join(t)

def good(p, D=None):
    s = str(p); digs = set(s)
    if D is not None and len(digs) != D: return False
    if not isprime(p): return False
    for a in digs:
        f = F(a, s)
        if f == s or not isprime(int(f)): return False
    return True

for line in open("writeup/data/p795_solutions.txt"):
    if line.startswith("#"): continue
    L, D, p = map(int, line.split())
    assert len(str(p)) == L and good(p, D), (L, D, p)
    assert p < 2**64 or pari.isprime(p) == 1
    print(f"L={L} D={D} p={p}: ok;", "transforms:", " ".join(f"F({a})={F(a, str(p))}" for a in sorted(set(str(p)))))

# Vantieghem's larger numbers (not proven minimal)
for p, D in ((102756857398639201, 9), (102344856798765493201, 10)):
    assert good(p, D) and pari.isprime(p) == 1
    for a in sorted(set(str(p))): assert pari.isprime(int(F(a, str(p)))) == 1
    print(f"Vantieghem D={D} p={p}: ok (p and all transforms proven prime with PARI)")

# lemma: a digit that occurs once is a fixed point of F; with 10 digits each twice the digit sum is 90, divisible by 3
assert F('7', '1234567') == '1234567'
assert sum(2 * d for d in range(10)) % 3 == 0
print("lemma ok: every used digit occurs at least twice; a pandigital solution needs at least 21 digits")

# brute force L = 7 (D = 3) and L = 8 (D = 4): smallest solution among all numbers with D distinct digits each >= 2 times
def brute(L, D):
    best = None
    for digs in itertools.combinations(range(10), D):
        # count vectors: each >= 2, total L
        for extra in itertools.combinations_with_replacement(range(D), L - 2 * D):
            ms = []
            for k, d in enumerate(digs): ms += [str(d)] * (2 + extra.count(k))
            for perm in multiset_permutations(ms):
                if perm[0] == '0' or perm[-1] not in '1379': continue
                v = int("".join(perm))
                if best is not None and v >= best: continue
                if good(v, D): best = v
    return best
assert brute(7, 3) == 1051051 and brute(8, 4) == 39703709
print("brute force L=7, 8 agrees with the ordered search")
print("all checks passed")
