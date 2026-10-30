"""Puzzle 1035: x = tau(sigma(x)) * sigma(tau(x)) only for x = 1, 468, 3240 (proof in results/1035.md)."""
import math
from fractions import Fraction
import numpy as np
from sympy import divisor_count, divisor_sigma, primerange

# 1. The constant: tau(m) <= c * m^(1/4) with c^4 = prod_{p < 16} max_a (a+1)^4 / p^a.
c4 = Fraction(1)
for p in primerange(2, 16):
    c4 *= max(Fraction((a + 1)**4, p**a) for a in range(60))
assert c4 == Fraction(127401984, 25025)
# primes p >= 17 contribute 1: (a+1)^4 <= 16^a <= p^a
assert all((a + 1)**4 <= 16**a for a in range(200))
c = float(c4) ** 0.25 * (1 + 1e-12)

# 2. The bound f(x) = c^2 sqrt(x) (1+ln x)^(1/4) (1 + ln c + ln(x)/4) is below x at x = 10^6,
#    and f(x)/x decreases (derivative of its log in t = ln x is < 1/4 + 1/4 - 1/2 = 0).
f = lambda x: c * c * math.sqrt(x) * (1 + math.log(x))**0.25 * (1 + math.log(c) + math.log(x) / 4)
assert f(1e6) / 1e6 < 0.93 and math.log(c) > 0

# 3. Sanity checks of the two lemmas on m <= 10^5.
for m in range(1, 100001):
    t, s = int(divisor_count(m)), int(divisor_sigma(m))
    assert t**4 <= c4 * m and s <= m * (1 + math.log(m)) + 1e-9

# 4. Exhaustive below 10^6, two ways.
X = 10**6
S = 8 * X  # sigma(x) < 8x below 10^6 (sigma(x)/x < 5 there)
tau = np.zeros(S + 1, dtype=np.int32); sig = np.zeros(X + 1, dtype=np.int64)
for d in range(1, S + 1):
    tau[d::d] += 1
for d in range(1, X + 1):
    sig[d::d] += d
assert sig[1:].max() <= S
sol_np = [x for x in range(1, X) if tau[sig[x]] * int(sig[tau[x]]) == x]
sol_sympy = [x for x in range(1, X) if divisor_count(divisor_sigma(x)) * divisor_sigma(divisor_count(x)) == x]
assert sol_np == sol_sympy == [1, 468, 3240], (sol_np, sol_sympy)
print("solutions below 10^6:", sol_np, "; bound ratio at 10^6:", round(f(1e6) / 1e6, 4))
print("all checks pass")
