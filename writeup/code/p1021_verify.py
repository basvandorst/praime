"""Puzzle 1021. Least prime p(k) such that
  A: p(k) + p(k+1) + 1          B: p(k) + p(k+1) + p(k+2)
is squarefree with exactly n prime factors, for n = 1..NMAX.

Method. The sum x is an increasing function of k, so the least k belongs to the least admissible x.
For n >= 3 every admissible x is odd and squarefree with omega(x) = n (p = 2 gives 6 and 10, two factors).
The odd squarefree numbers with n prime factors are generated in exact increasing order with a heap over
index tuples (i_1 < ... < i_n) into the odd primes; each x is tested against the primes nearest (x-1)/2
(type A) or x/3 (type B) (Volpatti's method).  The first x that works is the least, no bound is needed.
Primes proved with PARI isprime (APR-CL).  Output: data/p1021_table.txt.  Needs gmpy2, sympy, cypari2."""
import sys, time, heapq
from gmpy2 import mpz, next_prime, is_prime
from sympy import primerange, factorint
import cypari2
pari = cypari2.Pari(); pari.allocatemem(1 << 30)     # APR-CL on 250-digit numbers needs a large stack
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 44
primes = [mpz(p) for p in primerange(3, 200000)]

def prev_prime(n):
    n = mpz(n) - 1
    while not is_prime(n, 25): n -= 1
    return n
def check_A(x):
    p = prev_prime((x - 1) // 2 + 1); q = next_prime(p)
    return (p, q) if p + q + 1 == x else None
def check_B(x):
    p1 = prev_prime(x // 3 + 1); p0 = prev_prime(p1); p2 = next_prime(p1); p3 = next_prime(p2)
    for t in ((p0, p1, p2), (p1, p2, p3)):
        if sum(t) == x: return t
    return None
def squarefree_in_order(n):
    """yield odd squarefree numbers with exactly n prime factors in increasing order"""
    start = tuple(range(n)); x = mpz(1)
    for i in start: x *= primes[i]
    heap = [(x, start)]; seen = {start}
    while heap:
        x, t = heapq.heappop(heap)
        yield x
        for j in range(n):                                  # raise index j by one if it stays strictly increasing
            nj = t[j] + 1
            if j + 1 < n and nj >= t[j + 1]: continue
            u = t[:j] + (nj,) + t[j + 1:]
            if u in seen: continue
            seen.add(u); heapq.heappush(heap, (x // primes[t[j]] * primes[nj], u))
def proven(n): return pari.isprime(int(n)) == 1

page = {  # least p(k) from the page (Cleary / Volpatti), n = 1..20
 'A': [5, 2, 191, 1201, 20521, 390001, 2424907, 205307593, 4071444847, 74345698897,
       5444525719111, 87193673083987, 6428314689524191, 219136252305472951,
       15506131856082359923, 1052297628357522954853, 48549698841393869136463,
       2829235940486873030746081, 198422619164776743160946947, 19916760708314661900454948627],
 'B': [5, 2, 331, 3049, 35227, 299903, 7234517, 158055193, 1524467891, 146981019161,
       2924264157767, 124113516642097, 4535946133123387, 132476253972317237,
       7276195760630481613, 711930353432117784941, 37511148993323851080139,
       3130872412489286082022393, 185807294787963924506776321, 11071065345213582903890991271]}
assert check_A(mpz(13)) == (5, 7) and check_B(mpz(23)) == (5, 7, 11)       # n = 1: 13 and 23 prime
assert 2 + 3 + 1 == 6 == 2 * 3 and 2 + 3 + 5 == 10 == 2 * 5                 # n = 2 with p = 2
rows = [(1, 'A', 5, 13, '13'), (1, 'B', 5, 23, '23'), (2, 'A', 2, 6, '2 3'), (2, 'B', 2, 10, '2 5')]
t0 = time.time(); mism = []
for n in range(3, NMAX + 1):
    res = {}; cnt = 0
    for x in squarefree_in_order(n):
        cnt += 1
        for typ, chk in (('A', check_A), ('B', check_B)):
            if typ in res: continue
            r = chk(x)
            if r: res[typ] = (x, r, cnt)
        if len(res) == 2: break
    for typ in 'AB':
        x, pr, c = res[typ]
        f = factorint(int(x)); assert all(e == 1 for e in f.values()) and len(f) == n and all(proven(v) for v in pr)
        for u, v in zip(pr, pr[1:]): assert next_prime(u) == v        # consecutive primes
        if n <= 20 and pr[0] != page[typ][n - 1]: mism.append((typ, n, int(pr[0]), page[typ][n - 1]))
        rows.append((n, typ, int(pr[0]), int(x), " ".join(str(q) for q in sorted(f))))
    print(f"n={n} A: p={res['A'][1][0]} ({res['A'][2]} candidates)  B: p={res['B'][1][0]} ({res['B'][2]})  {time.time()-t0:.0f}s", flush=True)
    with open("writeup/data/p1021_table.txt", "w") as fh:          # rewritten after every n
        fh.write("# puzzle 1021: n type least_p sum prime_factors_of_sum\n")
        for r in rows: fh.write(f"{r[0]} {r[1]} {r[2]} {r[3]} {r[4]}\n")
print("mismatches with the page (type, n, ours, page):", mism)
print("table written for n <=", NMAX)
