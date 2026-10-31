# Puzzle 937: m(n) = X p_(n+1) - B, B = p_n#, X = least integer > B/p_(n+1) coprime to B.
# Claim: m(n) = least positive t*p_(n+1) - r (r = B mod p_(n+1)) with no prime factor <= p_n; so m(n) is 1 or prime
# whenever m(n) < p_(n+2)^2.  Checked for n <= 600 both ways (n <= 1500 takes 4 minutes, same result).
from math import gcd, log, exp
from sympy import prime, primorial, isprime
first = []
worst = 0
for n in range(1, 601):
    B = primorial(n); P1 = prime(n + 1); P2 = prime(n + 2)
    # definition
    X = B // P1 + 1
    while gcd(X, B) != 1: X += 1
    m = X * P1 - B
    # reformulation: smallest t >= 1 with t*P1 - r coprime to B
    r = B % P1; t = 1
    while gcd(t * P1 - r, B) != 1: t += 1
    assert m == t * P1 - r
    assert m == 1 or (isprime(m) and m > P1)          # 1 or a prime above p_(n+1)
    assert m < P2 * P2                                  # the condition that forces "1 or prime"
    worst = max(worst, t / (P2 * P2 / P1))
    if n <= 15: first.append(m)
assert first == [1, 19, 19, 43, 17, 43, 1, 31, 41, 43, 137, 199, 103, 59, 79]   # the page's values
print("largest t / (p_(n+2)^2/p_(n+1)) for n <= 600:", round(worst, 4))
print("all claims verified")
