"""Puzzle 455: checks for results/455.md."""
from sympy import primerange, divisor_count, divisor_sigma, primefactors, factorint

def holds(fac):
    """1 + tau(n) sigma(n) == n rad(n) for n = prod p^a, from the factorisation."""
    n = r = A = S = 1
    for p, a in fac:
        n *= p**a; r *= p; A *= a + 1; S *= (p**(a + 1) - 1) // (p - 1)
    return 1 + A * S == n * r

# 1. Prime powers: p^a works iff a = p - 2 (checked p < 200, a <= 300).
assert all(holds([(p, a)]) == (a == p - 2) for p in primerange(3, 200) for a in range(1, 301))
assert not any(holds([(2, a)]) for a in range(1, 300))

# 2. Two prime factors: the bounds of the proof leave finitely many (a, b); none works.
cases = 0
for a in range(1, 80):
    for b in range(1, 80):
        A = (a + 1) * (b + 1)
        if not (2**(a + 1) < 2 * A and 2**(b + 1) <= A * A + A):
            continue
        cases += 1
        ps = list(primerange(2, A + 1))
        for p in ps:
            for q in ps:
                if q != p and (p - 1) * (q - 1) < A:
                    assert not holds([(p, a), (q, b)]), (p, a, q, b)
assert max(a for a in range(1, 80) for b in range(1, 80)
           if 2**(a + 1) < 2 * (a + 1) * (b + 1) and 2**(b + 1) <= ((a + 1) * (b + 1))**2 + (a + 1) * (b + 1)) < 79
print("(a, b) cases with two prime factors:", cases)

# 3. Brute-force cross-check of the two-prime case on a box.
assert not any(holds([(p, a), (q, b)]) for p in primerange(2, 80) for q in primerange(2, 80) if p < q
               for a in range(1, 14) for b in range(1, 14))

# 4. Direct check with sympy for n <= 2*10^5, and the five numbers found by the C search to 10^10.
direct = [n for n in range(2, 200001)
          if 1 + divisor_count(n) * divisor_sigma(n) == n * __import__('math').prod(primefactors(n))]
assert direct == [3, 125, 16807], direct
for n in (3, 125, 16807, 2357947691):
    (p, a), = factorint(n).items()
    assert a == p - 2 and holds([(p, a)])
for n in (77, 128807153):      # n divides 1 + tau*sigma, but the equation fails
    t = 1 + divisor_count(n) * divisor_sigma(n)
    assert t % n == 0 and t != n * __import__('math').prod(primefactors(n))
print("all checks pass")
