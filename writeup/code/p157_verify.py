# Puzzle 157: there are exactly 20 perfect zip primes (writeup/157.md).
# Run: python3 p157_verify.py   (about 5 seconds; needs sympy)
from itertools import product
from sympy import isprime, factorint


def parts(s, k):
    """Zip the digit string s into k parts: part j gets the digits j, j+k, j+2k, ... (0-based)."""
    return [s[j::k] for j in range(k)]


def zippable(s, k):
    return all(isprime(int(p)) for p in parts(s, k))


def perfect(s):
    return all(zippable(s, k) for k in range(1, len(s) + 1))


def interleave(a, b):
    return ''.join(a[i // 2] if i % 2 == 0 else b[i // 2] for i in range(len(a) + len(b)))


# 0. The example on the page: 27239 is a zip prime of order 4 and not 5-zippable.
assert all(zippable('27239', k) for k in range(1, 5)) and not zippable('27239', 5)

# 1. The lemma: zipping A = s[0::2] into m parts gives the odd-numbered parts of s zipped into 2m parts,
#    and B = s[1::2] gives the even-numbered ones. It is pure index arithmetic; checked here on every
#    string of length <= 9 over the digits 2, 3, 5, 7.
for n in range(1, 10):
    for t in product('2357', repeat=n):
        s = ''.join(t)
        A, B = s[0::2], s[1::2]
        for m in range(1, n // 2 + 1):
            assert parts(A, m) == parts(s, 2 * m)[0::2] and parts(B, m) == parts(s, 2 * m)[1::2]

# 2. The recursion: every perfect zip prime of length n >= 2 interleaves perfect zip primes of lengths
#    ceil(n/2) and floor(n/2).  Build the candidates from the halves and keep the perfect ones.
P = {1: ['2', '3', '5', '7']}
cand = {}
for n in range(2, 15):
    cand[n] = [interleave(a, b) for a in P[(n + 1) // 2] for b in P[n // 2]]
    P[n] = [c for c in cand[n] if perfect(c)]
expected = {1: ['2', '3', '5', '7'], 2: ['23', '37', '53', '73'], 3: ['223', '233', '337', '523', '733', '773'],
            4: ['5233'], 5: ['33377', '72733'], 6: ['272333', '572333'], 7: ['5222333']}
for n in range(1, 8):
    assert P[n] == expected[n]

# 3. Lengths 8..14: 1, 2, 4, 4, 4, 2, 1 candidates (18 in all), none perfect.
assert [len(cand[n]) for n in range(8, 15)] == [1, 2, 4, 4, 4, 2, 1]
for n in range(8, 15):
    assert P[n] == []
    for c in cand[n]:
        k = next(k for k in range(1, n + 1) if not zippable(c, k))
        print('n = %2d  %s  fails at k = %d' % (n, c, k), dict(factorint(int(c))) if k == 1 else '(prime)')

# 4. Induction: for n >= 15 the half of length ceil(n/2) has length between 8 and n-1, so P[n] is empty.
#    Independent brute force over all numbers with digits 2, 3, 5, 7 (every digit of a perfect zip prime
#    is a 1-digit part, hence prime) up to 10 digits gives the same 20 numbers.
brute = [s for n in range(1, 11) for s in (''.join(t) for t in product('2357', repeat=n)) if perfect(s)]
assert brute == [x for n in range(1, 8) for x in P[n]]
print('perfect zip primes:', brute)
print('count:', len(brute), '- all checks passed')
