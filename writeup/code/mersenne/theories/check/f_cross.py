# cross-check F up to 1e5 with sympy.n_order
from sympy import primerange, n_order, isprime
c = 0
for q in primerange(3, 10**5 + 1):
    o = n_order(2, q)
    if o % 2 == 1 and isprime(o):
        c += 1
print("F <=1e5 via n_order:", c)
