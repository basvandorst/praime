"""Independent J_inf: the local quantity depends only on the equality pattern of the 4 points 0, x1=a, x2=a+b, x3=a+b+c mod p
(x1,x2,x3 uniform independent). Enumerate the 15 set partitions of {0,1,2,3}; count of realisations = (p-1)(p-2)...(p-k+1).
Cross-check by brute force over (a,b,c) mod p for small p, and check the marginal identities E_c f4 = f3, E_a f4 = f3(b,c)."""
import math, itertools
from fractions import Fraction
import numpy as np
import sympy

def partitions(s):
    if not s:
        yield []; return
    first, rest = s[0], s[1:]
    for p in partitions(rest):
        for i in range(len(p)):
            yield p[:i] + [[first] + p[i]] + p[i+1:]
        yield [[first]] + p

PARTS = list(partitions([0, 1, 2, 3]))
assert len(PARTS) == 15

def f(nu, k, p):
    return (1 - nu / p) / (1 - 1 / p) ** k

def J_p(p):
    tot = 0.0; wsum = 0.0
    for P in PARTS:
        k = len(P)
        if k > p: continue
        cnt = math.prod(p - i for i in range(1, k))
        blk = {e: i for i, B in enumerate(P) for e in B}
        n = lambda idx: len({blk[e] for e in idx})
        f4 = f(k, 4, p); f3ab = f(n([0, 1, 2]), 3, p); f3bc = f(n([1, 2, 3]), 3, p); f2 = f(n([1, 2]), 2, p)
        w = cnt / p ** 3
        wsum += w * f4
        if f4 > 0:
            tot += w * f4 * math.log(f4 * f2 / (f3ab * f3bc))
    return tot, wsum

def brute(p):
    tot = 0.0; marg_err = 0.0
    for a in range(p):
        for b in range(p):
            nb = len({0, b % p}); f2 = f(nb, 2, p)
            f3ab = f(len({0, a, (a + b) % p}), 3, p)
            s = 0.0
            for c in range(p):
                f4 = f(len({0, a, (a + b) % p, (a + b + c) % p}), 4, p)
                f3bc = f(len({0, b, (b + c) % p}), 3, p)
                s += f4
                if f4 > 0:
                    tot += f4 * math.log(f4 * f2 / (f3ab * f3bc))
            marg_err = max(marg_err, abs(s / p - f3ab))
    return tot / p ** 3, marg_err

for p in (3, 5, 7, 11, 13):
    jp, ws = J_p(p); jb, me = brute(p)
    print(p, 'partition', jp, 'brute', jb, 'mean f4', ws, 'marginal err', me)

# all odd primes up to 1e7 exactly by partition formula; scaling check
ps = list(sympy.primerange(3, 10 ** 6))
J = np.array([J_p(p)[0] for p in ps])
pa = np.array(ps, float)
for q in (101, 199, 401, 1009, 10007, 99991):
    i = ps.index(q); print('p', q, 'J_p', J[i], 'p^3 J_p', J[i] * q ** 3, 'p^4 J_p', J[i] * q ** 4)
head400 = J[pa <= 400].sum(); tot = J.sum()
# tail beyond 1e6 using p^-3 or p^-4 asymptote from last values
c3 = J[-1] * pa[-1] ** 3
print('sum p<=400', head400, 'sum p<1e6', tot, 'tail >400', tot - head400)
print('J_3 J_5 J_7', J[0], J[1], J[2])
print('J_inf nats ~', tot, 'bits', tot / math.log(2))
