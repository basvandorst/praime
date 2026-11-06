"""Puzzle 285: K(2^n - 1) mod n, where K is Kempner's function (least m with N | m!).

For every n <= NMAX (default 200) and every multiple of 8 up to 400: factor 2^n - 1 with PARI, compute K exactly
from the prime powers, take the largest prime factor p and d = ord_p(2), t = n/d, and check the law
    p = 1 (mod d);  p mod n = 1 + d*(((p-1)/d) mod t);  (p-1)/d odd  <=>  p = +-3 (mod 8);  d odd => 2d | p-1,
its corollaries for t = 1, 2, 4, and the page's observations (residue 1 except n in {28, 52, 68, 92, 124, 156, 172,
244, ...} where it is n/2 + 1; Andrew Rupinski's rule p mod 4 for n = 4 mod 8).
"""
import sys
import cypari2
from sympy import n_order, isprime

pari = cypari2.Pari()
pari.allocatemem(10**9)

def kempner_pp(q, e):              # least m with q^e | m!  (Legendre's formula, m a multiple of q)
    m = q
    while True:
        v, x = 0, m
        while x:
            x //= q; v += x
        if v >= e:
            return m
        m += q

def factor(N):
    f = pari.factor(N)
    return [(int(f[0][i]), int(f[1][i])) for i in range(len(f[0]))]

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 200
rows = []
for n in list(range(2, NMAX + 1)) + list(range(NMAX + 8 - NMAX % 8, 401, 8)):
    N = 2**n - 1
    fa = factor(N)
    K = max(kempner_pp(q, e) for q, e in fa)
    p, e = max(fa)
    assert K == p and e == 1                         # Kempner value = largest prime factor, which is simple
    d = n_order(2, p); t = n // d
    assert n % d == 0 and (p - 1) % d == 0
    r = p % n
    assert r == (1 + d * (((p - 1) // d) % t)) % n
    assert (((p - 1) // d) % 2 == 1) == (p % 8 in (3, 5))
    if d % 2 == 1:
        assert (p - 1) % (2 * d) == 0
    if t == 1:
        assert r == 1
    if t == 2:
        assert r == (1 if p % 8 in (1, 7) else n // 2 + 1)
    if t == 4 and d % 2 == 1:
        assert r == (1 if p % 4 == 1 else n // 2 + 1)
    if n % 8 == 4 and t in (1, 2, 4):
        assert r == (1 if p % 4 == 1 else n // 2 + 1)      # Rupinski's rule, proved for t <= 4
    if n % 4 == 2 and t in (1, 2):
        assert r == 1
    rows.append((n, p, d, t, r))

exc = [(n, r, t) for n, p, d, t, r in rows if r != 1]
print("residue != 1:", exc)
assert [n for n, r, t in exc if n <= NMAX] == [28, 52, 68, 92, 124, 156, 172][:len([n for n, r, t in exc if n <= NMAX])]
assert all(r == n // 2 + 1 and n % 8 == 4 for n, r, t in exc)
print("t = n/ord_p(2) by n mod 8:", {k: sorted({t for n, p, d, t, r in rows if n % 8 == k}) for k in range(8)})
print("multiples of 8 with t = 2:", [(n, p % 8, r) for n, p, d, t, r in rows if n % 8 == 0 and t == 2])
print(f"checked n <= {NMAX} and all multiples of 8 up to 400")
