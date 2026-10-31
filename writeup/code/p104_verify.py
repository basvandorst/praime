# Puzzle 104 (sum of the squares of the two halves): the number of solution pairs of length L is d(10^L + 1)/2 - 1.
# Run: python3 p104_verify.py   (about 10 seconds; needs sympy)
from math import isqrt
from sympy import factorint, divisor_count

# 1. De Geest's counts on the page (L = 2..18) against d(10^L+1)/2 - 1, and against Nash's 2^(F-1) - 1.
page = {2: 0, 4: 1, 6: 1, 8: 1, 10: 3, 12: 3, 14: 7, 16: 15, 18: 3}
for L in range(2, 31, 2):
    M = 10**L + 1
    f = factorint(M)
    assert all(p % 4 == 1 for p in f)                        # every prime factor of b^2 + 1 (b even) is 1 mod 4
    pairs = divisor_count(M) // 2 - 1
    if L in page:
        assert pairs == page[L]
    assert all(e == 1 for e in f.values()) and pairs == 2**(len(f) - 1) - 1   # squarefree for L <= 30
    print('L = %2d  pairs = %3d   10^L+1 = %s' % (L, pairs, ' * '.join(map(str, sorted(f)))))


# 2. Brute force: x b + y = x^2 + y^2 (0 < x < b, b = 10^m) iff 4x(b - x) + 1 = (2y - 1)^2.
def solutions(m):
    b = 10**m
    out = []
    for x in range(1, b):
        t = 4 * x * (b - x) + 1
        s = isqrt(t)
        if s * s == t:
            y = (s + 1) // 2
            assert x * b + y == x * x + y * y
            out.append((x, y))
    return out


for m in range(1, 7):
    sol = solutions(m)
    assert len(sol) == 2 * (divisor_count(10**(2 * m) + 1) // 2 - 1)
    assert all((10**m - x, y) in sol for x, y in sol)            # the pairing x <-> b - x, same y

# 3. The page's L = 18 list: 999999000_001000000 is a solution; its partner is 000001000_001000000,
#    not 000001000_999000000 as printed.
b = 10**9
assert 999999000 * b + 1000000 == 999999000**2 + 1000000**2
assert 1000 * b + 1000000 == 1000**2 + 1000000**2
assert 1000 * b + 999000000 != 1000**2 + 999000000**2

# 4. 10^L + 1 need not be squarefree: 17 || 10^8 + 1, so by lifting the exponent 17^2 || 10^136 + 1, and then
#    d(10^136+1)/2 - 1 = 3 * 2^(F-2) - 1 > 2^(F-1) - 1: Nash's "at most 2^(F-1) - 1" is reversed.
assert (10**8 + 1) % 17 == 0 and (10**8 + 1) % 17**2 != 0
assert (10**136 + 1) % 17**2 == 0 and (10**136 + 1) % 17**3 != 0
print('all checks passed')
