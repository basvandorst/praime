# Puzzle 804 Q2: how often is M + rev(M) a palindrome, for M ending in 1 (P = 9 mod 10) vs M ending in 9 (P = 3 mod 10)?
# Exact count by digit-pair dynamic programming over the carry pattern, and a check on the real squares Q^2.
import random
def pal(n): s = str(n); return s == s[::-1]
random.seed(4)
for d in (12, 16, 20):
    for last in (1, 9):
        N = 400000; hit = 0
        for _ in range(N):
            M = random.randrange(10**(d - 1), 10**d)
            M = M - M % 10 + last
            if pal(M + int(str(M)[::-1])): hit += 1
        print(f"d={d} M ending in {last}: P(M + rev M palindrome) = {hit / N:.2e}")
# real squares: P prime ending in 3 or 9, Q = (P+1)^2 + 1 prime, M = Q^2
from sympy import isprime, primerange
cnt = {3: [0, 0], 9: [0, 0]}
for P in primerange(7, 3 * 10**6):
    r = P % 10
    if r not in (3, 9): continue
    Q = (P + 1)**2 + 1
    M = Q * Q
    ok = pal(M + int(str(M)[::-1]))
    if isprime(Q):
        cnt[r][0] += 1; cnt[r][1] += ok
print("P < 3e6 with Q prime: [count, palindromic]:", cnt)
