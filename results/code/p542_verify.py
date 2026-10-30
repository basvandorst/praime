"""Puzzle 542: independent checks for results/542.md (plain sympy/gmpy2, no C)."""
import gmpy2
from sympy import divisor_sigma, primerange, isprime, perfect_power

def odd_prime_power(s):
    """(q, k) if s = q^k with q an odd prime, else None."""
    if s % 2 == 0:
        return None
    if isprime(s):
        return (s, 1)
    pp = perfect_power(s)
    if pp and isprime(pp[0]):
        return pp
    return None

# 1. Direct: sigma(n) is odd only for n = m^2 or 2m^2. Check all of them up to 2*10^10.
hits = []
for m in range(1, 100001):
    for n in (m*m, 2*m*m):
        r = odd_prime_power(int(divisor_sigma(n)))
        if r and r[1] > 1:
            hits.append(n)
assert sorted(hits) == [81, 400], hits

# 2. Lemma: sigma(2^a) = 2^(a+1) - 1 is never q^j with j >= 2 (checked a <= 200; proof in the write-up).
assert all(not (odd_prime_power(2**(a+1) - 1) or (0, 1))[1] > 1 for a in range(1, 201))

# 3. Q2 with q = 2: sigma(p^m) is never a power of 2 for m >= 2 (checked p < 1000, m <= 30).
assert not any(gmpy2.is_power(int(divisor_sigma(p**m))) and int(divisor_sigma(p**m)) & (int(divisor_sigma(p**m)) - 1) == 0
               for p in primerange(3, 1000) for m in range(2, 31))

# 4. Q2 beyond the 10^18 search for small p: sigma(p^m) = (p^(m+1)-1)/(p-1) a perfect power
#    only for 3^4 (= 11^2) and 7^3 (= 20^2), for p < 1000 and 2 <= m <= 400.
pw = [(p, m) for p in primerange(2, 1000) for m in range(2, 401)
      if gmpy2.is_power((p**(m+1) - 1)//(p - 1))]
assert pw == [(3, 4), (7, 3)], pw
print("all checks pass")
