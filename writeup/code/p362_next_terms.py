"""Puzzle 362: the greedy sequence of odd primes with prod (1 - 1/a(n)) = 1/2.

a(1) = 3; with r(n) = prod_{i<=n} (1 - 1/a(i)) the next term is the least prime
q with r(n) (1 - 1/q) > 1/2, i.e. q > r(n) / (r(n) - 1/2).  This is OEIS A339344.
Prints a(1..12) exactly; each term is proved prime with PARI's isprime (APR-CL/ECPP).
"""
from fractions import Fraction
import gmpy2, cypari2, sys
pari = cypari2.Pari()
pari.allocatemem(2 * 10**9)
K = int(sys.argv[1]) if len(sys.argv) > 1 else 12
r = Fraction(1)
seq = []
for n in range(1, K + 1):
    if n == 1:
        q = 3
    else:
        bound = r / (r - Fraction(1, 2))          # need q > bound
        q = int(gmpy2.next_prime(int(bound)))     # least prime > floor(bound) = least prime > bound (bound is not an integer unless prime>bound anyway)
        assert q > bound
        # check no prime between bound and q
        assert int(gmpy2.next_prime(int(bound))) == q
        assert r * (1 - Fraction(1, q)) > Fraction(1, 2)
        # minimality: the prime before q does not work
        pq = int(gmpy2.prev_prime(q)) if hasattr(gmpy2, 'prev_prime') else int(pari.precprime(q - 1))
        assert r * (1 - Fraction(1, pq)) <= Fraction(1, 2)
    proved = bool(pari.isprime(q))
    assert proved
    r *= (1 - Fraction(1, q))
    seq.append(q)
    d = len(str(q))
    print(f"a({n}) = {q if d <= 70 else str(q)[:30] + '...' + str(q)[-30:]}  ({d} digits, proved prime: {proved})")
    print(f"   excess r(n) - 1/2 = {float(r - Fraction(1,2)):.3e}")
with open(sys.argv[2] if len(sys.argv) > 2 else '/dev/null', 'w') as f:
    for n, q in enumerate(seq, 1):
        f.write(f"{n} {q}\n")
