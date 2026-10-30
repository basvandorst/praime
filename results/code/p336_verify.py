"""Puzzle 336: sigma(n) < prime(n) for every n > 4 (proof in results/336.md)."""
import math
from fractions import Fraction
from sympy import divisor_sigma, prime, primerange

# 1. The crossover: n(1 + ln n) <= n(ln n + ln ln n - 1)  iff  ln ln n >= 2  iff  n >= e^(e^2) = 1618.17...
N0 = math.exp(math.exp(2))
assert 1618 < N0 < 1619 and math.log(math.log(1619)) >= 2

# 2. Lemma: sigma(n)/n <= H_n < 1 + ln n for n >= 2 (checked exactly for n <= 2000 with fractions for H_n).
H = Fraction(0)
for n in range(1, 2001):
    H += Fraction(1, n)
    if n >= 2:
        assert Fraction(int(divisor_sigma(n)), n) <= H and float(H) < 1 + math.log(n)

# 3. Direct check up to 10^5 (covers 5 <= n <= 1618 with room to spare), and the exceptions n = 2, 4.
P = list(primerange(2, 1300000))
assert len(P) >= 100000
fails = [n for n in range(1, 100001) if divisor_sigma(n) >= P[n - 1]]
assert fails == [2, 4], fails
assert divisor_sigma(2) == prime(2) == 3 and divisor_sigma(4) == prime(4) == 7

# 4. Dusart's bound p_n > n(ln n + ln ln n - 1) holds on our range (sanity check of the cited theorem).
assert all(P[n - 1] > n * (math.log(n) + math.log(math.log(n)) - 1) for n in range(2, 100001))
print("all checks pass")
