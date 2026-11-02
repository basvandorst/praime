# Puzzle 1281 (writeup/1281.md): checks the primes on the page and ours. For a prime with blocks of k digits:
# every block is a k-digit number with no zero, the blocks are pairwise coprime, and (cases 3, 4) no block is prime.
# Also recomputes the block-count bounds: a composite block below 10^k has a prime factor below 10^(k/2), and all the
# blocks' prime factors are distinct.
# Run: python3 p1281_verify.py   (a second; needs sympy)
from math import gcd, isqrt
from sympy import isprime, primerange, factorint

def check(name, n, k, noprimes):
    s = str(n); assert isprime(n), name
    blocks = [int(s[i:i + k]) for i in range(len(s) - k + 1)]
    assert all(b >= 10 ** (k - 1) and '0' not in str(b) for b in blocks), name
    assert all(gcd(a, b) == 1 for i, a in enumerate(blocks) for b in blocks[i + 1:]), name
    if noprimes: assert not any(isprime(b) for b in blocks), name
    print('%-10s %-46s prime, %d digits, %d blocks pairwise coprime%s' % (name, s, len(s), len(blocks), ', none prime' if noprimes else ''))

check('case 1', 8997941735371319, 2, False)
check('case 2', 88779773974337937335997199113713191933173139, 3, False)
check('case 3', 94951, 2, True)
check('case 4', 97377998951, 3, True)
for k in (2, 3):
    small = list(primerange(2, isqrt(10 ** k - 1) + 1))
    assert all(min(factorint(b)) in small for b in range(10 ** (k - 1), 10 ** k) if not isprime(b))
    print('k = %d: every composite block has a prime factor among %s, so at most %d composite blocks, %d digits' % (k, small, len(small), len(small) + k - 1))
print('all checks passed')
