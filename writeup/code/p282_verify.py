"""Puzzle 282: check every claim from scratch with plain integers.
Needs sympy and cypari2 (PARI's isprime(n, 2) is an APR-CL primality proof)."""
import random
import cypari2
from sympy import isprime, totient, divisor_sigma

pari = cypari2.Pari()
proven = lambda n: bool(pari.isprime(n, 2))
rev = lambda n: int(str(n)[::-1])

PHI = [7198858264383196858264206778858264384978840444384977076264384978957]
SIGMA = [98639,
         986146894593861385310540531145946885405310540613861468945938639,
         9861468854061468945946885406138613861385311459468945938531145938639,
         199437576265342819943782825448418985342819943782825448418957626282591]

for p in PHI + SIGMA:
    r = rev(p)
    assert p != r and proven(p)
    m = 19 if r % 19 == 0 else 47
    b = r // m
    assert r == m * b and proven(b) and b % m
    if p in PHI:
        assert totient(p) == p - 1 == (m - 1) * (b - 1) == totient(r)
        assert m * rev(r) == totient(m) * r - m * (totient(m) - 1)
    else:
        assert divisor_sigma(p) == p + 1 == (m + 1) * (b + 1) == divisor_sigma(r)
        assert m * rev(r) == divisor_sigma(m) * r + m * (divisor_sigma(m) - 1)
    print(len(str(p)), 'digits, rev(p) =', m, '* prime:', 'phi' if p in PHI else 'sigma', 'equal')

for _ in range(10000):
    n = random.randrange(1, 10**random.randrange(2, 40))
    if n % 10 == 0:
        continue
    k = len(str(n))
    assert (rev(n) - n) % 9 == 0
    assert (rev(n) - (-1) ** (k - 1) * n) % 11 == 0
print('rev(n) = n (mod 9) and rev(n) = (-1)^(k-1) n (mod 11) on 10000 random n')
