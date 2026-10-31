# Puzzle 41 F1: A a 9-digit anagram of 1..9, |A - rev(A)| palindromic. Reduce to d_i = a_(8-i) - a_i, i = 0..3.
from itertools import permutations, product
from collections import defaultdict
pal = lambda n: str(n) == str(n)[::-1]
# (1) brute force over all 9! anagrams
vals = defaultdict(int)
for t in permutations('123456789'):
    A = int(''.join(t)); D = abs(A - int(''.join(t[::-1])))
    if pal(D): vals[D] += 1
print("all anagrams: palindromic |A - rev A| =", dict(vals))
# (2) the reduction: |D| = 99 |1010101 d0 + 101010 d1 + 10100 d2 + 1000 d3|
for t in permutations('123456789', 9):
    a = [int(c) for c in t[::-1]]; A = int(''.join(t))
    d = [a[8 - i] - a[i] for i in range(4)]
    assert abs(A - int(''.join(t[::-1]))) == 99 * abs(1010101 * d[0] + 101010 * d[1] + 10100 * d[2] + 1000 * d[3])
    break
# (3) which difference patterns give palindromes (all d in [-8,8]\{0}), and are they realisable with distinct digits?
def realisable(d):
    for pairs in product(*[[(x, x + di) for x in range(1, 10) if 1 <= x + di <= 9] for di in d]):
        used = [v for p in pairs for v in p]
        if len(set(used)) == 8: return True
    return False
pats = defaultdict(list)
for d in product([x for x in range(-8, 9) if x], repeat=4):
    D = 99 * abs(1010101 * d[0] + 101010 * d[1] + 10100 * d[2] + 1000 * d[3])
    if pal(D) and D > 0: pats[D].append(d)
for D, L in sorted(pats.items()):
    ok = [d for d in L if realisable(d)]
    print(f"{D}: {len(L)} difference patterns, {len(ok)} realisable, e.g. {L[:4]}")
# (4) question I: A/9 palindromic -> always 8 digits, hence divisible by 11
from sympy import isprime
q9 = [int(''.join(t)) // 9 for t in permutations('123456789') if pal(int(''.join(t)) // 9)]
assert len(q9) == 34 and all(len(str(q)) == 8 and q % 11 == 0 for q in q9)
print("I: 34 palindromes A/9, all with 8 digits, all divisible by 11")
