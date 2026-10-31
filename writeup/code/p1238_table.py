from sympy import primerange, n_order
from p1238_criterion import revdigits, simulate
import json
rows=[]; NMAX=120
for p in primerange(7, 1000):
    T = n_order(10, p); H = {pow(10,k,p) for k in range(T)}
    alg = any((1-u) % p in H for u in H)          # some u, v in <10> with u + v = 1
    L = NMAX + 3*T + 10; b, Y = revdigits(p, L); best=None
    for N in range(2, NMAX+1):
        for J in range(1, N):
            if (10**N + 10**J - 1) % p == 0 and (10**N + 10**J - 1)//p + 10**J*(Y % 10**(N-J)) < 10**N:
                assert simulate(b, N, J, L) == b
                best=(N,J); break
        if best: break
    rows.append((p,T,alg,best, b[:best[0]] if best else None))
json.dump(rows, open('rev_table.json','w'))
n=len(rows); print(n, "primes; algebraic solutions:", sum(r[2] for r in rows), "; valid reverse rule with N<=%d:"%NMAX, sum(1 for r in rows if r[3]))
print("algebraic but no valid rule found:", [r[0] for r in rows if r[2] and not r[3]][:30])
for r in rows[:16]: print(r)
