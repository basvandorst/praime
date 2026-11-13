"""Second, separate check of the shorter prime-representing constants (exact rational arithmetic, sympy primes).
Builds each constant from its definition as an exact fraction, truncates it to d decimal digits, runs the
recurrence with exact fractions and counts the primes it returns before the first wrong one."""
import math
from fractions import Fraction as Fr
from sympy import prime, nextprime

def radix(kind, p):
    if kind == 'fridman': return p
    if kind == 'new': return math.floor(math.log(p) ** 2) + 4
    if kind == 'half': return math.floor(math.log(p) ** 2 / 2) + 3

def build(kind, terms=1500):
    start = 3 if kind == 'half' else 2
    f = Fr(0); den = 1; p = start
    for _ in range(terms):
        q = nextprime(p); g = q - p
        R = radix(kind, p)
        assert (g // 2 if kind == 'half' else g) + 1 <= R
        den *= R; f += Fr(g // 2 if kind == 'half' else g, den); p = q
    return (1 if kind == 'half' else 2) + f

def decode(kind, x, limit=2000):
    out = [2] if kind == 'half' else []
    p_true = 3 if kind == 'half' else 2
    for _ in range(limit):
        m = math.floor(x)
        p = 2 * m + 1 if kind == 'half' else m
        if p != p_true: break
        out.append(p)
        x = m + radix(kind, p) * (x - m)
        p_true = nextprime(p_true)
    return out

for kind in ('fridman', 'new', 'half'):
    c = build(kind)
    res = []
    for d in (20, 100, 1000):
        x = Fr(math.floor(c * 10 ** (d - 1)), 10 ** (d - 1))
        got = decode(kind, x)
        res.append(f'{d} digits -> {len(got)} primes')
    print(f'{kind}: ' + '; '.join(res))
