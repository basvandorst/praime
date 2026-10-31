# Puzzle 759: at most 10 palindromic primes in arithmetic progression (writeup/759.md).
# Run: python3 p759_verify.py   (about 20 seconds; needs sympy, and cypari2 for the primality proofs of section 1)
from sympy import isprime, primerange

pal = lambda n: str(n) == str(n)[::-1]

# 1. The 1999 record: 10 palindromic primes in AP, (p, d) = (742950290870000078092059247, 1010100000000000).
p, d = 742950290870000078092059247, 1010100000000000
terms = [p + k * d for k in range(10)]
assert all(pal(t) for t in terms) and not pal(p + 10 * d) and not pal(p - d)
try:
    import cypari2
    pari = cypari2.Pari()
    assert all(int(pari('isprime(%d, 1)' % t)) == 1 for t in terms)       # proof (Pocklington-Lehmer / APRCL)
    print('1999 record: 10 proven palindromic primes')
except ImportError:
    assert all(isprime(t) for t in terms)
    print('1999 record: 10 palindromic (BPSW) probable primes (install cypari2 for proofs)')
# Its structure is the extremal case of step 4: lowest nonzero digit of d at position j = 11 with e = 1,
# and digit 11 of the terms runs through 0, 1, ..., 9; all terms have 27 digits and the same first digit.
j = len(str(d)) - len(str(d).rstrip('0'))
assert j == 11 and str(d)[len(str(d)) - 1 - j] == '1'
assert [t // 10**j % 10 for t in terms] == list(range(10))
assert len({len(str(t)) for t in terms}) == 1 and len({str(t)[0] for t in terms}) == 1

# 2. The page's examples (Seiji, Jan, Emmanuel) are palindromic primes in AP, except Emmanuel's first one:
#    (7190917, 10101) gives 7201018; with d = 101010 (Jan's pair) it works.
examples = [(900010009, 10202010, 5), (7190917, 101010, 5), (10872427801, 1000100010, 5),
            (1012717172101, 101010101010, 5), (101270555072101, 10101000101010, 5), (92134043129, 1000100010, 6),
            (1156021206511, 101010000, 6), (7272836382727, 101010101010, 6), (105433737334501, 10101000101010, 6),
            (118023131320811, 100010001000, 7), (349263515362943, 1010100000, 7), (15519112621191551, 10101010101000, 7),
            (13133112721133131, 101010001010100, 7), (13032102720123031, 101010001010100, 8),
            (1701910118110191071, 11011110111101100, 8), (1863101311009001131013681, 11011110111101100000, 9),
            (19125452191, 10101000, 6), (16313711111731361, 1000100010000, 7), (169979300121003979961, 1010100000000, 8)]
for a, b, n in examples:
    assert all(pal(a + k * b) and isprime(a + k * b) for k in range(n)), (a, b, n)
    assert b % 10 == 0 or n < 6
assert not pal(7190917 + 10101)

# 3. Exhaustive check of the whole statement for small numbers: among all palindromic primes below 10^7,
#    the longest AP has length <= 10; every AP of length >= 6 has 10 | d and constant length.
P = [q for q in primerange(2, 10**7) if pal(q)]
S = set(P)
longest = 0
for i, a in enumerate(P):
    for b in P[i + 1:]:
        dd = b - a
        n = 2
        while a + n * dd in S:
            n += 1
        if n >= 6:
            assert dd % 10 == 0 and len({len(str(a + k * dd)) for k in range(n)}) == 1
        longest = max(longest, n)
print('palindromic primes below 10^7:', len(P), ' longest AP among them:', longest)
assert longest <= 10

# 4. The key step for general sizes is index arithmetic: with 10^j || d (lowest nonzero digit e of d at position j),
#    floor(P_k / 10^j) = floor(p / 10^j) + k d / 10^j, so digit j of P_k is (p_j + k e) mod 10 and the digits
#    below j never change. Checked here on random large APs.
import random
rnd = random.Random(1)
for _ in range(2000):
    jj = rnd.randrange(1, 15)
    ee = rnd.randrange(1, 10)
    dd = (rnd.randrange(0, 10**6) * 10 + ee) * 10**jj
    a = rnd.randrange(10**20, 10**21)
    for k in range(12):
        t = a + k * dd
        assert t % 10**jj == a % 10**jj and t // 10**jj % 10 == (a // 10**jj % 10 + k * ee) % 10
print('all checks passed')
