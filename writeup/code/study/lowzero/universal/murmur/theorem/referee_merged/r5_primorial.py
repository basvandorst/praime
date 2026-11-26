# Referee: class-number means at a primorial-like level N = 17# = 510510 (omega = 7), own code.
import math, sys
import numpy as np
from sympy import primerange
import cypari2
pari = cypari2.Pari(); pari.allocatemem(2 * 10**8)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 510510
from sympy import factorint
f = factorint(N); assert all(e == 1 for e in f.values()); PN = sorted(f)
P = list(primerange(2, 10 ** 6))
def u(l): return 1.0 / ((l - 1) ** 2 * (l + 1))
def e(l, r): return l * l / (l * l - 1.0) if r % l == 0 else 1 - u(l)
def E(r, excl):
    v = 1.0
    for l in P:
        if l not in excl: v *= e(l, r)
    return v
ps = list(primerange(N + 1, int(1.05 * N)))
print("N = %d = %s, %d primes in (N, 1.05N]" % (N, PN, len(ps)))
for r in [0, 1, 2, 3]:
    vals = []
    for p in ps:
        D = N * (4 * p - r * r * N)
        if D <= 0: continue
        vals.append(math.pi * float(pari.qfbhclassno(D)) / math.sqrt(D))
    vals = np.array(vals)
    mu = vals.mean(); se = vals.std() / math.sqrt(len(vals))
    own = E(r, set(PN)); prime_level = E(r, set())
    print(" r=%d (%d p): mean L* = %.4f (se %.4f) | own local prediction %.4f (z %+.2f) | prime-level nu'(r) %.4f (z %+.1f)"
          % (r, len(vals), mu, se, own, (mu - own) / se, prime_level, (mu - prime_level) / se))
