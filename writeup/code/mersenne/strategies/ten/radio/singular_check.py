"""Check of the identity 2 C2 S(n) = sum_q mu(q)^2 c_q(n) / phi(q)^2 for even n (c_q = Ramanujan sum), partial sums to q <= 2*10^5."""
from sympy import totient, mobius, primefactors
from math import gcd
C2 = 0.6601618158468696
def c(q, n):
    d = q // gcd(q, n); return int(mobius(d)) * int(totient(q)) // int(totient(d))
for n in (4, 6, 10, 30, 210, 1000):
    s = sum(c(q, n) / int(totient(q)) ** 2 for q in range(1, 200001) if mobius(q) != 0)
    S = 1.0
    for p in primefactors(n):
        if p > 2: S *= (p - 1) / (p - 2)
    print('n = %4d: series %.5f, 2 C2 S(n) = %.5f' % (n, s, 2 * C2 * S))
