"""Share of primes p <= X (p not dividing a) for which a is a primitive root, against Artin's constant with the
Heilbronn-Hooley correction for square-free a != -1:  A if a != 1 mod 4,  A (1 - mu(|a|) prod_{q | a} 1/(q^2 - q - 1))
if a = 1 mod 4.  A = prod_q (1 - 1/(q(q-1))).  Predictions written before the run: base 5 lands near 0.3936, base -3 near
0.4488, base 21 near 0.3721, bases 2, 3, 6, 7, 10 near 0.3740; base 4 (a square) only at p = 3.   usage: artin.py X out.txt"""
import sys, math
import numpy as np
from sympy import factorint, mobius, primerange
X = int(sys.argv[1]); out = open(sys.argv[2], 'w')
def pr(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')
spf = np.zeros(X + 1, np.int64)
for i in range(2, X + 1):
    if spf[i] == 0: spf[i::i][spf[i::i] == 0] = i
P = [int(p) for p in np.nonzero(spf == np.arange(X + 1))[0] if p >= 2]
def pf(n):
    s = set()
    while n > 1: q = int(spf[n]); s.add(q); n //= q
    return s
A = 1.0
for q in primerange(2, 10 ** 7): A *= 1 - 1 / (q * (q - 1))
pr('Artin constant from the primes below 10^7: %.10f (literature 0.3739558136)' % A)
F = {p: pf(p - 1) for p in P}
for a in (2, 3, 5, 6, 7, 10, 13, 21, -3, 4):
    n = c = 0
    for p in P:
        if a % p == 0: continue
        n += 1; g = a % p
        if all(pow(g, (p - 1) // q, p) != 1 for q in F[p]): c += 1
    if a == 4: pred = 0.0
    elif a % 4 == 1 or (a < 0 and (a % 4) == 1):
        m = int(mobius(abs(a))); corr = 1.0
        for q in factorint(abs(a)): corr *= 1 / (q * q - q - 1)
        pred = A * (1 - m * corr)
    else: pred = A
    sd = math.sqrt(pred * (1 - pred) / n) if 0 < pred < 1 else 0
    pr('a = %3d: %7d of %7d primes, share %.5f; predicted %.5f; difference %+.5f (%.1f binomial sd)' % (a, c, n, c / n, pred, c / n - pred, (c / n - pred) / sd if sd else float('nan')))
