"""Puzzle 870: check the three nested prime magic squares 4x4 of grade 2 in results/870.md."""
from sympy import isprime

def is_prime_magic(sq):
    n = len(sq)
    lines = [list(r) for r in sq] + [[sq[i][j] for i in range(n)] for j in range(n)]
    lines += [[sq[i][i] for i in range(n)], [sq[i][n - 1 - i] for i in range(n)]]
    flat = [x for r in sq for x in r]
    return len({sum(l) for l in lines}) == 1 and len(set(flat)) == n * n and all(isprime(x) for x in flat)

def nest(sq, L, R):
    """Puzzle 869 rule m = (L*10^LN + n)*10^LR + R; L or R may be absent (None)."""
    LN = max(len(str(x)) for r in sq for x in r)
    f = lambda n: int((str(L) if L is not None else '') + str(n).zfill(LN) + (str(R) if R is not None else ''))
    out = [[f(x) for x in r] for r in sq]
    # same as the arithmetic form of the rule
    LR = len(str(R)) if R is not None else 0
    assert all(f(x) == ((L or 0) * 10**LN + x) * 10**LR + (R or 0) for r in sq for x in r)
    return out

cases = [  # (inner square, L, R)
    ([[191, 239, 367, 541], [499, 409, 197, 233], [379, 331, 401, 227], [269, 359, 373, 337]], None, 3),
    ([[163, 557, 631, 863], [733, 761, 271, 449], [701, 463, 719, 331], [617, 433, 593, 571]], 3, None),
    ([[139, 281, 577, 859], [727, 709, 191, 229], [593, 409, 547, 307], [397, 457, 541, 461]], 2, 1),
]
for inner, L, R in cases:
    outer = nest(inner, L, R)
    assert is_prime_magic(inner) and is_prime_magic(outer)
    print(sum(inner[0]), sum(outer[0]), outer)
print("all checks pass")
